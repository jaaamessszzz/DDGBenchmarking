#Removes unnecessary entries in json file based on existing PDB_REDO pdb's in a directory
import json
import os

def main():
    #Initialize stuff
    jsoninfo  = json.load(open("/kortemmelab/home/james.lucas/160429-kyleb_zemu-jl-brub-rscr-unminjump/data/blank_job_dict.json"))
    predIDs = os.listdir("/kortemmelab/home/james.lucas/Zemu-PDB_REDO_Dataset/data")

    json_set = set()
    PDB_REDO_set = set()

    for json_predID in jsoninfo:
        json_set.add(json_predID)

    for predID in predIDs:
        PDB_REDO_set.add(predID)

    difference = json_set - PDB_REDO_set

    #Deletes unwanted entries from jsoninfo
    for item in difference:
        jsoninfo.pop(item)

    #Write updated jsoninfo to file
    open("blank_job_dict_updated.json", "w").write(
        json.dumps(jsoninfo, sort_keys=True,separators=(',', ': ')))

main()