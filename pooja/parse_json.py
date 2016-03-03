##
import json
from pprint import pprint
import sys
from itertools import izip
import os
import numpy
import random

with open('ZEMu.export.json') as data_file:
    zemu = json.load(data_file)

zemu_list=[]

for i in range(len(zemu['Data'])):
    zemu_resid=str(zemu['Data'][i]['PDBMutations'][0]['ResidueID']).strip()

    zemu_list.append(str(zemu['Data'][i]['PDBMutations'][0]['ID'])+" "+str(zemu['Data'][i]['PDBMutations'][0]['PDBFileID'])+" "+str(zemu['Data'][i]['PDBMutations'][0]['Chain'])+ str(zemu['Data'][i]['PDBMutations'][0]['WildTypeAA'])+ zemu_resid+ str(zemu['Data'][i]['PDBMutations'][0]['MutantAA'])+" "+str(zemu['Data'][i]['PublishedDDG']))

#generate random numbers and keep in list
with open("rand_ints") as f:
    num=[int(line.split(" ")[0]) for line in f]

#print a random entry from the zemu dataset using the random number list
for i in range(50,70):
    zemu_id=str(zemu_list[num[i]]).split(" ")[0]
    pdb=str(zemu_list[num[i]]).split(" ")[1]
    mut=str(zemu_list[num[i]]).split(" ")[2]
    pub=str(zemu_list[num[i]]).split(" ")[3]
    print zemu_id, pdb, mut, pub

#copy pdb from database, copy contents, python resfile_find_neighbors.py
for i in range(50,70):
    pdb=str(zemu_list[num[i]]).split(" ")[1]
    mut=str(zemu_list[num[i]]).split(" ")[2]
    if not os.path.exists(pdb):
        os.mkdir(pdb)
    os.chdir("%s/"%pdb)
    if not os.path.exists(mut):
        os.mkdir(mut)
    #os.mkdir("%s"%mut)
    os.chdir("%s"%mut)
    os.system("cp /netapp/database/pdb/remediated/pdb/%s/pdb%s.ent.gz ." %(pdb[1:3].lower(), pdb.lower()))
    os.system("gunzip *.gz")
    os.system("mv *.ent %s.pdb"%pdb)
    os.system("cp ~/DDG/new_runs/contents/* .")
    print pdb, mut
    os.system("python resfile_find_neighbors.py %s.pdb %s 8"%(pdb, mut))
    os.chdir("../../")

#submit_multistate.py (this script submits backrub and repack jobs)
for i in range(50, 70):
    pdb=str(zemu_list[num[i]]).split(" ")[1]
    mut=str(zemu_list[num[i]]).split(" ")[2]
    try:
        os.chdir("%s/%s"%(pdb,mut))
        os.system("cp ~/DDG/new_runs/contents/* .")
        os.system("python submit_multistate.py %s.pdb %s"%(pdb, mut))
        os.chdir("../../")
    except:
        pass

#submit_final_min.py (this script submits the last minimization step for every member in the ensemble, post backrub and repack)
for i in range(0,30):
    pdb=str(zemu_list[num[i]]).split(" ")[0]
    mut=str(zemu_list[num[i]]).split(" ")[1]
    try:
        os.chdir("%s/%s"%(pdb,mut))
        os.system("cp ~/DDG/new_runs/contents/* .")
        os.system("python submit_final_min.py")
        os.chdir("../../")
    except:
        pass

#separate monomers and partners, and then rescore using fixbb
for i in range(0,50):
    pdb=str(zemu_list[num[i]]).split(" ")[0]
    mut=str(zemu_list[num[i]]).split(" ")[1]
    try:
        os.chdir("%s/%s"%(pdb,mut))
        os.system("cp ~/DDG/new_runs/contents/* .")
        os.system("python submit_separate_partners.py")
        os.system("python submit_final_natro.py")
        os.chdir("../../")
    except:
        pass

#compute DDGs affinities
for i in range(37,40):
    pdb=str(zemu_list[num[i]]).split(" ")[0]
    mut=str(zemu_list[num[i]]).split(" ")[1]
    pub=str(zemu_list[num[i]]).split(" ")[2]
    print pdb, mut, pub
    try:
        os.chdir("%s/%s/"%(pdb,mut))
        os.system("cp ~/DDG/new_runs/contents/* .")
        os.chdir("wildtype/")
        os.system("python find_ddg_af_wt.py %s.pdb %s %s"%(pdb,mut,pub))
    except:
        pass
