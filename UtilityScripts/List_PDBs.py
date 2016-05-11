#Prints out all PDB files in subdirectories of workingdir
#List_PDBs.py
import os
import csv

workingdir = '/Users/jameslucas/Kortemme_Rotation/PDB_REDO'
#workingdir = '/netapp/home/james.lucas/160315-kyleb_james-backrub-rscript/data/'
#workingdir = '/Users/jameslucas/Kortemme_Rotation/output/'

pdbcsv = open('PDB_list.csv', 'wb')
csvwriter = csv.writer(pdbcsv, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)

pdb_set = set()
for i in os.listdir(workingdir):
    if i.endswith('.pdb'):
        pdb_set.add(i)
    else:
        continue

for pdbid in pdb_set:
    print pdbid
    csvwriter.writerow([pdbid])