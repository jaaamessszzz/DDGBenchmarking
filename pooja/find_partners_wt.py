import os
import sys

#####
#separate monomer and partner after final minimization
#####

mut_complexes=[]
wt_complexes=[]

mut=os.getcwd()[38:43]
#pdb=os.getcwd()[33:37]
#print mut
for filename in os.listdir(os.getcwd()):
    if filename.endswith("0001_0001_last_mutate_0001_0001.pdb"):
        mut_complexes.append(filename)
    if filename.endswith("0001_0001_last_repack_0001_0001.pdb"):
        wt_complexes.append(filename)

mut_complexes.sort()
wt_complexes.sort()

#pdb_complex=open(sys.argv[1], "r")
#PDB= sys.argv[1].split(".")[0]
#monomer=open("%s_monomer.pdb"%PDB, "w")
#partner=open("%s_partner.pdb"%PDB, "w")

#json_mutation=str(sys.argv[1])
chain_id=str(mut[0])
#print chain_id

for wt_complex in wt_complexes:
    open_wt=open(wt_complex, "r")
    monomer=open("%s_monomer.pdb"%str(wt_complex).rstrip('.pdb'), "w")
    partner=open("%s_partner.pdb"%str(wt_complex).rstrip('.pdb'), "w")
    for line in open_wt:
        if "ATOM" in line[0:4]:
            if chain_id in line[21:22]:
                monomer.write(str(line))
            if chain_id not in line[21:22]:
                partner.write(str(line))

for mut_complex in mut_complexes:
    open_mut=open(mut_complex, "r")
    monomer=open("%s_monomer.pdb"%str(mut_complex).rstrip('.pdb'), "w")
    partner=open("%s_partner.pdb"%str(mut_complex).rstrip('.pdb'), "w")
    for line in open_mut:
        if "ATOM" in line[0:4]:
            if chain_id in line[21:22]:
                monomer.write(str(line))
            if chain_id not in line[21:22]:
                partner.write(str(line))
