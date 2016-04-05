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

        for line in enumerate(open(filename, 'r')):
            #Finds first instance of 'fa_atr' and adds values to fattydict line-by-line until it reaches the dashed cutoff line
            if line[1].find("fa_atr") == 1:
                if counter % 2 == 0:
                    struct_type = 'WT'
                else:
                    struct_type = 'Mutant'

                linecounter = 0
                currentline = line[0] + 1

                if structID not in fattydict[i]['structIDs']:
                    fattydict[i]['structIDs'][structID] = {}
                fattydict[i]['structIDs'][structID][struct_type] = {}

                try:
                    while linecache.getline(filename, currentline).strip() != '-----------------------------------------':
                        scores = linecache.getline(filename, currentline)
                        parsed_scores = scores.split()
                        fattydict[i]['structIDs'][structID][struct_type][parsed_scores[0]] = float( parsed_scores[1] )
                        linecounter = linecounter + 1
                        currentline = currentline + 1
                    sumscore = linecache.getline(filename, line[0] + linecounter + 2)
                    parsed_sumscore = sumscore.split()
                    fattydict[i]['structIDs'][structID][struct_type][parsed_sumscore[0]] = float( parsed_sumscore[2] )
                    counter = counter + 1
                except:
                    print "Oops, something went wrong here..."

            elif line[1].find("reported success in") > 1:
                timeline = line[1].split()
                fattydict[i]['structIDs'][structID]['Runtime'] = float( timeline[5] ) / 60.0 # Convert seconds to minutes
                structID = structID +1
            else:
                continue

        if verbose:
            print str(i) + ": " + str(structID - 1) + " structures completed"

        #Parse output file for Max VMem usage (GB), start/end times, and return code
        files = os.listdir( os.path.join(workingdir, str(i) ) )
        for doc in files:
            if doc.startswith("SubmitRun_DDG_Zemu_General.py.o"):
                output = doc
            else:
                continue

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
                fattydict[i]['Max virtual memory usage'] = float(mem_usage)

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
