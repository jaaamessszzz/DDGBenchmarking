#Prints out all PDB files in subdirectories of cwd
#List_PDBs.py
import os
import csv

cwd = os.getcwd()
print cwd

pdbcsv = open('PDB_list.csv', 'wb')
csvwriter = csv.writer(pdbcsv, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)

pdb_set = set()
for dirpath, dirnames, filenames in os.walk(cwd):
    for filename in filenames:
        if filename.endswith('.pdb'):
            pdb_set.add(filename)
            
print pdb_set

for pdbid in pdb_set:
    csvwriter.writerow([pdbid])
