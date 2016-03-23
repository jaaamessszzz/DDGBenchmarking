#Iterates through directories to open Rosetta.out and creates Nested Dictionaries:
#>Prediction ID
#>>Structure ID (nstruct 1-100)
#>>>Structure Type (WT, Mut, DDG)
#>>>>Score Type

import os
import linecache

def parse_rosetta_out(workingdir):
    fattydict = {}
    for i in os.listdir(workingdir):
        if os.path.isdir(workingdir + i):
            fattydict[i] = {}
            structID = 1
            counter = 0
            filename = workingdir + i + "/rosetta.out"
            for line in enumerate(open(filename, 'r')):
                if line[1].find("fa_atr") == 1:
                    if counter % 2 == 0:
                        linecounter = 0
                        currentline = line[0] + 1
                        fattydict[i]['WT_' + str(structID)] = {}
                        while linecounter < 15:
                            scores = linecache.getline(filename, currentline)
                            parsed_scores = scores.split()
                            fattydict[i]['WT_' + str(structID)][parsed_scores[0]] = parsed_scores[1]
                            linecounter = linecounter + 1
                            currentline = currentline + 1
                        sumscore = linecache.getline(filename, line[0] + 17)
                        parsed_sumscore = sumscore.split()
                        fattydict[i]['WT_' + str(structID)][parsed_sumscore[0]] = parsed_sumscore[2]
                        counter = counter + 1
                    else:
                        linecounter = 0
                        currentline = line[0] + 1
                        fattydict[i]['Mutant_' + str(structID)] = {}
                        while linecounter < 15:
                            scores = linecache.getline(filename, currentline)
                            parsed_scores = scores.split()
                            fattydict[i]['Mutant_' + str(structID)][parsed_scores[0]] = parsed_scores[1]
                            linecounter = linecounter + 1
                            currentline = currentline + 1
                        sumscore = linecache.getline(filename, line[0] + 17)
                        parsed_sumscore = sumscore.split()
                        fattydict[i]['Mutant_' + str(structID)][parsed_sumscore[0]] = parsed_sumscore[2]
                        counter = counter + 1
                        structID = structID +1
                else:
                    continue
            print str(i) + ": " + str(structID - 1) + " structures completed"
    return fattydict

my_working_directory = '/netapp/home/james.lucas/DDG_Zemu_Backrub_Output'
parsed_dict = parse_rosetta_out(my_working_directory)
os.chdir(my_working_directory)

open("DDG_Data.json", "w").write(
    json.dumps(parsed_dict, sort_keys=True,separators=(',', ': ')))