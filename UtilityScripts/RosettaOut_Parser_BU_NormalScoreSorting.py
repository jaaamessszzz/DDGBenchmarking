#Iterates through directories to open Rosetta.out and creates Nested Dictionaries:
#>Prediction ID
#>>Structure ID (nstruct 1-100)
#>>>Structure Type (WT, Mut, DDG)
#>>>>Score Type

import os
import linecache
import json
from klab.Reporter import Reporter
import datetime

def parse_rosetta_out(workingdir, verbose = True):
    fattydict = {}
    unfinished = {}

    task_dirs = [os.path.join(workingdir, d) for d in os.listdir(workingdir) if os.path.isdir( os.path.join(workingdir, d) )]
    r = Reporter('parsing task directories', entries = 'directories')
    r.set_total_count( len(task_dirs) )

    #For each subdirectory (PredictionID) in the working directory (output/)...
    for task_dir in task_dirs:
        try:
            i = long(os.path.basename(task_dir))
        except ValueError:
            print 'Directory %s is not a number' % i
            continue
        
        fattydict[i] = {}
        fattydict[i]['structIDs'] = {}
        structID = 1
        counter = 0
        filename = os.path.join(task_dir, "rosetta.out")
        
        if not os.path.isfile( filename ):
            print 'Missing output file:', filename
            continue
        
        temp_dict = {}
        temp_dict['WT'] = {}
        temp_dict['Mutant'] = {}
        score_list = ["fa_atr",
                      "fa_rep",
                      "fa_sol",
                      "fa_intra_rep",
                      "mm_bend",
                      "fa_elec",
                      "hbond_sr_bb",
                      "hbond_lr_bb",
                      "hbond_bb_sc",
                      "hbond_sc",
                      "rama",
                      "omega",
                      "fa_dun",
                      "p_aa_pp",
                      "yhh_planarity",
                      "ref"]
        
        print i
        
        for line in enumerate(open(filename, 'r')):
            #WT or Mutant scores for current structID
            if counter % 2 == 0:
                struct_type = 'WT'
            else:
                struct_type = 'Mutant'
            
            if line[1].split()[0] == "Unbound":
                Total_scores = False
            if line[1].split()[0] == "Bound":
                Total_scores = False
            if line[1].split()[0] == "Scores":
                Total_scores = True
                
            #Looks for scoretype as first phrase in line, adds score to dict if present
            if Total_scores == True:
                for score in score_list:
                    if score in line[1].split()[0]:
                        parsed_scores = line[1].split()
                        temp_dict[struct_type][parsed_scores[0]] = float( parsed_scores[1] )

                if "Sum ddg: " in line[1]:
                    parsed_sumscore = line[1].split()
                    temp_dict[struct_type][parsed_sumscore[0]] = float( parsed_sumscore[2] )
                    counter = counter + 1
                    if counter % 2 == 0:
                        if len(temp_dict['Mutant']) != 0:
                            fattydict[i]['structIDs'][structID] = temp_dict
                            temp_dict = {}
                            temp_dict['Mutant'] = {}
                            temp_dict['WT'] = {}

                if "reported success in" in line[1]:
                    timeline = line[1].split()
                    fattydict[i]['structIDs'][structID]['Runtime'] = float( timeline[5] )
                    structID = structID + 1

        if verbose:
            print str(i) + ": " + str(structID - 1) + " structures completed"

        #Parse output file for Max VMem usage (GB), start/end times, and return code
        files = os.listdir(os.path.join(workingdir,  str(i)))
        for doc in files:
            if ".py.o" in doc:
                output = doc

                date_format_string = '%Y-%m-%d %H:%M:%S'
                outputfile = os.path.join(workingdir, os.path.join(str(i), output))
                for line in enumerate(open(outputfile, 'r')):
                    if line[1].startswith("Starting time:"):
                        starting_time = line[1][15:].strip()
                        fattydict[i]['Starting time'] = datetime.datetime.strptime(starting_time, date_format_string)
                    elif line[1].startswith("Ending time:"):
                        ending_time = line[1][13:].strip()
                        fattydict[i]['Ending time'] = datetime.datetime.strptime(ending_time, date_format_string)
                    elif line[1].startswith("Task return code:"):
                        task_return_code = line[1].split()
                        fattydict[i]['Task return code'] = long(task_return_code[3])
                    elif line[1].startswith("Max virtual memory usage:"):
                        mem_usage_parsed = line[1].split()
                        mem_usage = float(mem_usage_parsed[4][:-1])
                        size = mem_usage_parsed[4][-1:]
                        if size == 'M':
                            mem_usage = mem_usage / 1000
                        elif size == 'G':
                            continue
                        else:
                            print "Memory usage not measured in MB or GB!"
                        fattydict[i]['Max virtual memory usage:'] = float(mem_usage)

        #Keeps track of unfinished jobs
        if structID - 1 < 50: ###Change for each run!!!!
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
    print parsed_dict                
    #os.chdir(my_working_directory)

    #open("DDG_Data.json", "w").write(json.dumps(parsed_dict, sort_keys=True,separators=(',', ': ')))

if __name__ == '__main__':
    main()
