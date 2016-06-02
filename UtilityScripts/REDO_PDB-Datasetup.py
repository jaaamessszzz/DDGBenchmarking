#Copies subdirectories into new data directory, ignores .pdb's
#Copies and replaces REDO_PDB .pdb's into new data directory
import os
import shutil

#jsoninfo = json.load(open("/netapp/home/james.lucas/PDB_REDO_Stripped/blank_job_dict_updated.json"))
source = '/kortemmelab/home/james.lucas/160429-kyleb_zemu-jl-brub-rscr-unminjump/data/'
dest = '/kortemmelab/home/james.lucas/Zemu-PDB_REDO_Dataset/data/'
PDB_REDO_dir = '/kortemmelab/home/james.lucas/DDGBenchmarks_Test/PDB_REDO_Stripped/'

PDB_REDO_list = []

for PDB_REDO_files in os.listdir(PDB_REDO_dir):
    PDB_REDO_list.append(PDB_REDO_files)

counter = 0

for predID in os.listdir(source):
    if os.path.isdir(os.path.join(source, predID)):
        for data_item in os.listdir(os.path.join(source, predID)):
            for list_entry in PDB_REDO_list:
                if data_item == list_entry:
                    counter = counter +1
                    print os.path.join(dest, predID, list_entry)
                    #Copy Resfiles and stuff from source to PDB_REDO data directory
                    shutil.copytree(source + predID, dest + predID, ignore = shutil.ignore_patterns('*.pdb'))
                    #Copy PDB_REDO pdb's from PDB_REDO directory to new data directory
                    shutil.copy2(os.path.join(PDB_REDO_dir, list_entry), os.path.join(dest, predID, list_entry))

print counter

# for i in jsoninfo:
#     parsed = re.sub("/", ' ', jsoninfo[i]["input_file_list"][0])
#     data, filenum, pdbfile = parsed.split()
#     #Copy Resfiles and stuff from source to PDB_REDO data directory
#     shutil.copytree(source + i, dest + i, ignore = shutil.ignore_patterns('*.pdb'))
#     #Copy PDB_REDO pdb's from PDB_REDO directory to new data directory
#     shutil.copy2(PDB_REDO + "/" + pdbfile, dest + i + "/" + pdbfile)