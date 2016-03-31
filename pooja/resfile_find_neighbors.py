from Bio.PDB import *
import sys
import os
import shutil
import pprint

#####
#This script does 4 things. 1. It find neighbors (within an 8 angstrom shell of the mutation); 2. It generates resfiles and adds a list of neighbors to it; 3. It submits the initial minimization; and 4. It converts pdb numbers to rosetta numbers for the list of neighbor residues and saves it to a file called pivot_res.txt (The list written to this file will later be used for backrub)
#####

open_id = 'Open'

open_filename = str(sys.argv[1])
PDBID=open_filename.split(".")
#print PDBID[0]
parser = PDBParser(PERMISSIVE=1)

open_strct = parser.get_structure(open_id, open_filename)

chain_list=[]
neighbors= []
pdb_num = []
for model in open_strct:
    chain_number=0
    for chain in model:
        #print chain
        chain_number=chain_number+1
    for num in range(0,chain_number):
        chain_list.append(chr(num + ord('A')))

json_mutation=str(sys.argv[2])
##############
chain_id=str(json_mutation[0])
res_id=str(json_mutation[2:-1])
wt_aa=str(json_mutation[1])
mut_aa=str(json_mutation[-1])
##############

mut_chain = str(chain_id)
mut_pos = int(res_id)

#print mut_chain, mut_pos
mut_residue = open_strct[0][mut_chain][mut_pos]
#print mut_residue
for chain in chain_list:
    for residue in range(len(open_strct[0][chain])):
        residue += 1
        try:
            dist = mut_residue['CB'] - open_strct[0][chain][residue]['CB']
            if dist < float(sys.argv[3]):
                neighbors.append(str(residue)+" "+str(chain)+" "+"NATAA") #WANT CHAIN AS WELL
                if len(str(residue))==2:
                    pdb_num.append(str(chain)+"  "+str(residue)+" ")
                if len(str(residue))==3:
                    pdb_num.append(str(chain)+" "+str(residue)+" ")
                if len(str(residue))==4:
                    pdb_num.append(str(chain)+" "+str(residue))

                if len(str(residue+1))==2:
                    pdb_num.append(str(chain)+"  "+str(residue+1)+" ")
                if len(str(residue+1))==3:
                    pdb_num.append(str(chain)+" "+str(residue+1)+" ")
                if len(str(residue+1))==4:
                    pdb_num.append(str(chain)+" "+str(residue+1))

                if len(str(residue-1))==2:
                    pdb_num.append(str(chain)+"  "+str(residue-1)+" ")
                if len(str(residue-1))==3:
                    pdb_num.append(str(chain)+" "+str(residue-1)+" ")
                if len(str(residue-1))==4:
                    pdb_num.append(str(chain)+" "+str(residue-1))

        except KeyError:
            pass

#print neighbors
#print pdb_num

unique_pdb_num=[]
for elem in pdb_num:
    if elem not in unique_pdb_num:
        unique_pdb_num.append(elem)
#print unique_pdb_num

mutate_res=open("mutate.res", "w")
repack_res=open("repack.res", "w")
revert_res=open("revert.res", "w")

mutate_res.write("NATRO\nstart\n%s %s PIKAA %s\n" %(res_id, chain_id, mut_aa))
repack_res.write("NATRO\nstart\n")
revert_res.write("NATRO\nstart\n%s %s PIKAA %s\n" %(res_id, chain_id, wt_aa))

for i in neighbors:
    repack_res.write(str(i)+"\n")
    if i.split()[0] != res_id:
        mutate_res.write(str(i)+"\n")
        revert_res.write(str(i)+"\n")
mutate_res.close()
repack_res.close()
revert_res.close()
try:
    os.mkdir("wildtype")
    os.mkdir("mutation")
except:
    pass
os.system("cp repack.res revert.res mutation/")
os.system("cp repack.res mutate.res wildtype/")

mutation_suffix="_"+str(json_mutation[1:])
os.system("qsub submit_minimize.py %s %s"%(PDBID[0], mutation_suffix))

sys.path.insert(0, '/netapp/home/psuresh2/DDG')

from tools.fs.fsio import read_file
from tools.rosetta.map_pdb_residues import get_pdb_contents_to_pose_residue_map

pdb_content = read_file('./%s'%sys.argv[1])
rosetta_scripts_binary = '/netapp/home/thompson/Rosetta/rosetta_2014.22/source/bin/rosetta_scripts.linuxgccrelease'
success, result = get_pdb_contents_to_pose_residue_map(pdb_content, rosetta_scripts_binary, pdb_id = None, extra_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res')
#pprint.pprint(result)

rosetta_num=[]
for i in range(len(unique_pdb_num)):
    #print result[str(unique_pdb_num[i])]['pose_residue_id']
    rosetta_num.append(result[str(unique_pdb_num[i])]['pose_residue_id'])

#print rosetta_num

pivot_residues=open("pivot_res.txt", "w")
res=str(" ".join([str(x) for x in rosetta_num]))
pivot_residues.write(res)

minimize_pdb=open("pdblist.txt", "w")
minimize_pdb.write("%s_0001.pdb\n%s%s_0001.pdb"%(str(PDBID[0]), str(PDBID[0]), str(mutation_suffix)))
