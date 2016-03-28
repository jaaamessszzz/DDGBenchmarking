#Copies subdirectories into new data directory, ignores .pdb's
#Copies and replaces REDO_PDB .pdb's into new data directory
import os
import shutil
import json
import re

jsoninfo = json.load(open("/netapp/home/james.lucas/PDB_REDO/blank_job_dict_updated.json"))
source = '/netapp/home/james.lucas/160322-james-backrub-rscript-full/data/'
dest = '/netapp/home/james.lucas/Zemu-PDB_REDO_Dataset/data/'
PDB_REDO = '/netapp/home/james.lucas/PDB_REDO/'
for i in jsoninfo:
    parsed = re.sub("/", ' ', jsoninfo[i]["input_file_list"][0])
    data, filenum, pdbfile = parsed.split()
    #Copy Resfiles and stuff from source to PDB_REDO data directory
    shutil.copytree(source + i, dest + i, ignore = shutil.ignore_patterns('*.pdb'))
    #Copy PDB_REDO pdb's from PDB_REDO directory to new data directory
    shutil.copy2(PDB_REDO + "/" + pdbfile, dest + i + "/" + pdbfile)