#Iterates through directories to open Rosetta.out and creates Nested Dictionaries:
#>Prediction ID
#>>Structure ID (nstruct 1-100)
#>>>Structure Type (WT, Mut, DDG)
#>>>>Score Type

import os
import linecache
import json
from klab.Reporter import Reporter

def parse_rosetta_out(workingdir, verbose = True):
    fattydict = {}
    unfinished = {}

    task_dirs = [os.path.join(workingdir, d) for d in os.listdir(workingdir) if os.path.isdir( os.path.join(workingdir, d) )]
    r = Reporter('parsing task directories', entries = 'directories')
    r.set_total_count( len(task_dirs) )

    for task_dir in task_dirs:
        i = os.path.basename(task_dir)
        fattydict[i] = {}
        structID = 1
        counter = 0
        filename = os.path.join(task_dir, "rosetta.out")
        if not os.path.isfile( filename ):
            print 'Missing output file:', filename
            continue
        for line in enumerate(open(filename, 'r')):
            if line[1].find("fa_atr") == 1:
                if counter % 2 == 0:
                    linecounter = 0
                    currentline = line[0] + 1
                    fattydict[i]['WT_' + str(structID)] = {}
                    try:
                        while linecache.getline(filename, currentline).strip() != '-----------------------------------------':
                            scores = linecache.getline(filename, currentline)
                            parsed_scores = scores.split()
                            fattydict[i]['WT_' + str(structID)][parsed_scores[0]] = parsed_scores[1]
                            linecounter = linecounter + 1
                            currentline = currentline + 1
                        sumscore = linecache.getline(filename, line[0] + linecounter + 2)
                        parsed_sumscore = sumscore.split()
                        fattydict[i]['WT_' + str(structID)][parsed_sumscore[0]] = parsed_sumscore[2]
                        counter = counter + 1
                    except:
                        print "Oops, something went wrong here..."
                else:
                    linecounter = 0
                    currentline = line[0] + 1
                    fattydict[i]['Mutant_' + str(structID)] = {}
                    try:
                        while linecache.getline(filename, currentline).strip() != '-----------------------------------------':
                            scores = linecache.getline(filename, currentline)
                            parsed_scores = scores.split()
                            fattydict[i]['Mutant_' + str(structID)][parsed_scores[0]] = parsed_scores[1]
                            linecounter = linecounter + 1
                            currentline = currentline + 1
                        sumscore = linecache.getline(filename, line[0] + linecounter + 2)
                        parsed_sumscore = sumscore.split()
                        fattydict[i]['Mutant_' + str(structID)][parsed_sumscore[0]] = parsed_sumscore[2]
                        counter = counter + 1
                    except:
                        print "Oops, something went wrong here..."
            elif line[1].find("reported success in") > 1:
                timeline = line[1].split()
                fattydict[i]['Runtime_' + str(structID)] = timeline[5]
                structID = structID +1
            else:
                continue
        if verbose:
            print str(i) + ": " + str(structID - 1) + " structures completed"

        #Keeps track of unfinished jobs
        if structID - 1 < 100:
            unfinished[i] = structID - 1
        else:
            continue
        r.increment_report()
    r.done()
                
    return fattydict, unfinished

def main():
    my_working_directory = str(os.getcwd() + '/')
    print my_working_directory
    parsed_dict, unfinished_jobs = parse_rosetta_out(my_working_directory)

    print unfinished_jobs
    os.chdir(my_working_directory)

    open("DDG_Data.json", "w").write(json.dumps(parsed_dict, sort_keys=True,separators=(',', ': ')))

if __name__ == '__main__':
    main()
