# This script should be folded into the retirement phase of a job.
# It extracts the reported maxvmem from an output file and stores it in the database.

import sys
import os
import zipfile
import re
import subprocess
import time
import datetime

sys.path.insert(0, "../..")
sys.path.insert(0, "..")
from tools import colortext
from tools.fs.io import read_file, write_file
from tools.bio.pdb import PDB
from ddglib import ddgdbapi

def extract_data():
    ddGdb = ddgdbapi.ddGDatabase()
    #jobs = [r for r in ddGdb.execute_select('SELECT ID, maxvmem, DDGTime FROM Prediction WHERE Status="done" ORDER BY ID DESC')]
    jobs = [r for r in ddGdb.execute_select('SELECT ID, maxvmem, DDGTime FROM Prediction WHERE Status="done" ORDER BY ID DESC')]
    num_ids = len(jobs)
    count = 1
    results_root = '/kortemmelab/shared/DDG/jobs'
    for job in jobs:
        id = job['ID']
        maxvmem = job['maxvmem']
        DDGTime = job['DDGTime']
        if (maxvmem == 0 or maxvmem == None) or (DDGTime == 0 or DDGTime == None):
            colortext.message('%d (%d/%d)' % (id, count, num_ids))

            zipfile_path = os.path.join(results_root, '%d.zip' % id)
            try:
                z = zipfile.ZipFile(zipfile_path, 'r')
            except:
                colortext.error('MISSING FILE FOR %d' % PredictionID)
                continue
            file_list = z.namelist()

            stdout_files = [f for f in file_list if f.find('.cmd.o') != -1 and f.find('ddG') != -1]
            if not len(stdout_files) == 1:
                print('\n'.join(sorted(file_list)))
                assert(False)
            stdout = z.open(stdout_files[0], 'r').read()

            if stdout.find('<maxvmem>') != -1:
                maxvmem_in_GBs = []
                maxvmems = stdout[stdout.find('<maxvmem>') + 9 : stdout.find('</maxvmem>')].strip().split('\n')
                for maxvmem in maxvmems:
                    print('maxvmem: %s' % maxvmem)
                    mtchs = re.match('(^\d+[.]\d+)M$', maxvmem)
                    if mtchs:
                        maxvmem_in_GBs.append(float(mtchs.group(1))/1024)
                    else:
                        mtchs = re.match('(^\d+[.]\d+)G$', maxvmem)
                        if mtchs:
                            maxvmem_in_GBs.append(float(mtchs.group(1)))
                        else:
                            raise Exception('Could not parse "%s".' % maxvmem)
                if min(maxvmem_in_GBs)/max(maxvmem_in_GBs) >= 0.4:
                    maxvmem_in_GB = max(maxvmem_in_GBs)
                    ddGdb.execute('UPDATE Prediction SET maxvmem=%s WHERE ID=%s', parameters=(maxvmem_in_GB, id,))
            else:
                ddGdb.execute('UPDATE Prediction SET maxvmem=-1.0 WHERE ID=%s', parameters=(id,))

            if stdout.find('<startdate>') != -1 and stdout.find('<enddate>') != -1:
                startdate = stdout[stdout.find('<startdate>') + 11: stdout.find('</startdate>')].strip()
                enddate = stdout[stdout.find('<enddate>') + 9: stdout.find('</enddate>')].strip()
                startdate = datetime.datetime.strptime(startdate, '%a %b %d %H:%M:%S %Z %Y')
                enddate = datetime.datetime.strptime(enddate, '%a %b %d %H:%M:%S %Z %Y')
                DDGTime = float((enddate - startdate).total_seconds())/60.0
                ddGdb.execute('UPDATE Prediction SET DDGTime=%s WHERE ID=%s', parameters=(DDGTime, id,))
                print('DDG time: %0.2f' % DDGTime)

        count += 1

pdb_chain_lengths = {'1A23': {'A': 189},
 '1A43': {'A': 72},
 '1A5E': {'A': 156},
 '1AAR': {'A': 76, 'B': 76},
 '1ACB': {'E': 241, 'I': 63},
 '1AG2': {'A': 103},
 '1AJ3': {'A': 98},
 '1AKK': {'A': 104},
 '1AM7': {'A': 154, 'B': 154, 'C': 154},
 '1AMQ': {'A': 396},
 '1ANK': {'A': 214, 'B': 214},
 '1AON': {'A': 524,
          'B': 524,
          'C': 524,
          'D': 524,
          'E': 524,
          'F': 524,
          'G': 524,
          'H': 524,
          'I': 524,
          'J': 524,
          'K': 524,
          'L': 524,
          'M': 524,
          'N': 524,
          'O': 97,
          'P': 97,
          'Q': 97,
          'R': 97,
          'S': 97,
          'T': 97,
          'U': 97},
 '1APS': {'A': 98},
 '1AQH': {'A': 448},
 '1ARR': {'A': 53, 'B': 53},
 '1AV1': {'A': 201, 'B': 201, 'C': 201, 'D': 201},
 '1AXB': {'A': 263},
 '1AYE': {'A': 401},
 '1AYF': {'A': 103, 'B': 104},
 '1AZP': {'A': 66, 'B': 8, 'C': 8},
 '1B0O': {'A': 161},
 '1B26': {'A': 409, 'B': 409, 'C': 409, 'D': 409, 'E': 409, 'F': 409},
 '1B5M': {'A': 84},
 '1BCX': {'A': 185},
 '1BLC': {'A': 257},
 '1BNI': {'A': 108, 'B': 108, 'C': 107},
 '1BNL': {'A': 178, 'B': 178, 'C': 178, 'D': 178},
 '1BNZ': {'A': 64, 'B': 8, 'C': 8},
 '1BOY': {'A': 211},
 '1BP2': {'A': 123},
 '1BPI': {'A': 58},
 '1BRF': {'A': 53},
 '1BTA': {'A': 89},
 '1BVC': {'A': 153},
 '1C2R': {'A': 116, 'B': 116},
 '1C5G': {'A': 379},
 '1C9O': {'A': 66, 'B': 66},
 '1CAH': {'A': 258},
 '1CDC': {'A': 96, 'B': 96},
 '1CEY': {'A': 128},
 '1CHK': {'A': 238, 'B': 238},
 '1CLW': {'A': 543},
 '1CSP': {'A': 67},
 '1CYC': {'A': 103, 'B': 103},
 '1CYO': {'A': 88},
 '1DIL': {'A': 381},
 '1DIV': {'A': 149},
 '1DKT': {'A': 72, 'B': 71},
 '1DPM': {'A': 329, 'B': 329},
 '1EL1': {'A': 130, 'B': 130},
 '1FC1': {'A': 207, 'B': 207},
 '1FEP': {'A': 680},
 '1FKJ': {'A': 107},
 '1FLV': {'A': 168},
 '1FMK': {'A': 438},
 '1FNF': {'A': 368},
 '1FRD': {'A': 98},
 '1FTG': {'A': 168},
 '1FXA': {'A': 98, 'B': 98},
 '1G6N': {'A': 200, 'B': 201},
 '1H7M': {'A': 99},
 '1HFY': {'A': 120, 'B': 120},
 '1HFZ': {'A': 123, 'B': 122, 'C': 122, 'D': 124},
 '1HK0': {'X': 173},
 '1HME': {'A': 77},
 '1HNG': {'A': 175, 'B': 175},
 '1HTI': {'A': 248, 'B': 248},
 '1HUE': {'A': 90, 'B': 90},
 '1HZ6': {'A': 67, 'B': 63, 'C': 63},
 '1I5T': {'A': 104},
 '1IDS': {'A': 198, 'B': 198, 'C': 198, 'D': 198},
 '1IFB': {'A': 131},
 '1IGV': {'A': 75},
 '1IHB': {'A': 156, 'B': 156},
 '1IMQ': {'A': 86},
 '1IO2': {'A': 213},
 '1IOB': {'A': 153},
 '1IR3': {'A': 303, 'B': 6},
 '1IRL': {'A': 133},
 '1IRO': {'A': 53},
 '1JIW': {'I': 105, 'P': 470},
 '1K9Q': {'A': 40, 'B': 6},
 '1KDX': {'A': 81, 'B': 28},
 '1KFW': {'A': 435},
 '1L63': {'A': 162},
 '1LMB': {'1': 20, '2': 20, '3': 87, '4': 92},
 '1LS4': {'A': 164},
 '1LVE': {'A': 114},
 '1LZ1': {'A': 130},
 '1M7T': {'A': 107},
 '1MBG': {'A': 53},
 '1MGR': {'A': 97},
 '1MJC': {'A': 69},
 '1MSI': {'A': 66},
 '1MYL': {'A': 45, 'B': 40, 'C': 45, 'D': 45, 'E': 41, 'F': 40},
 '1N0J': {'A': 198, 'B': 198},
 '1OH0': {'A': 125, 'B': 128},
 '1OIA': {'A': 90, 'B': 86},
 '1ONC': {'A': 104},
 '1OTR': {'A': 49, 'B': 76},
 '1P2P': {'A': 124},
 '1PGA': {'A': 56},
 '1PHP': {'A': 394},
 '1PIN': {'A': 153},
 '1POH': {'A': 85},
 '1PQN': {'A': 127},
 '1QJP': {'A': 137},
 '1QLP': {'A': 372},
 '1QQV': {'A': 67},
 '1REX': {'A': 130},
 '1RGG': {'A': 96, 'B': 96},
 '1RHG': {'A': 145, 'B': 144, 'C': 145},
 '1RIS': {'A': 97},
 '1RN1': {'A': 103, 'B': 103, 'C': 104},
 '1ROP': {'A': 56},
 '1RRO': {'A': 108},
 '1RTB': {'A': 124},
 '1RTP': {'1': 109, '2': 109, '3': 109},
 '1RX4': {'A': 159},
 '1SAK': {'A': 42, 'B': 42, 'C': 42, 'D': 42},
 '1SCE': {'A': 97, 'B': 97, 'C': 111, 'D': 103},
 '1SHF': {'A': 59, 'B': 59},
 '1SHG': {'A': 57},
 '1SSO': {'A': 62},
 '1STN': {'A': 136},
 '1SUP': {'A': 275},
 '1TEN': {'A': 90},
 '1THQ': {'A': 147},
 '1TIT': {'A': 89},
 '1TPK': {'A': 88, 'B': 88, 'C': 88},
 '1TTG': {'A': 94},
 '1TUP': {'A': 196, 'B': 194, 'C': 195, 'E': 21, 'F': 21},
 '1U5P': {'A': 211},
 '1UBQ': {'A': 76},
 '1UWO': {'A': 91, 'B': 91},
 '1UZC': {'A': 69},
 '1VQB': {'A': 86},
 '1W4E': {'A': 45},
 '1W4H': {'A': 45},
 '1W99': {'A': 558},
 '1WIT': {'A': 93},
 '1WQ5': {'A': 258, 'B': 261},
 '1Y4Y': {'A': 87, 'B': 87, 'C': 87},
 '1YCC': {'A': 108},
 '1YEA': {'A': 112},
 '1ZNJ': {'A': 21,
          'B': 29,
          'C': 21,
          'D': 29,
          'E': 21,
          'F': 28,
          'G': 21,
          'H': 29,
          'I': 21,
          'J': 30,
          'K': 21,
          'L': 30},
 '219L': {'A': 164},
 '2A01': {'A': 243, 'B': 243, 'C': 243},
 '2A36': {'A': 59},
 '2ABD': {'A': 86},
 '2AC0': {'A': 199,
          'B': 198,
          'C': 195,
          'D': 198,
          'E': 12,
          'F': 12,
          'G': 12,
          'H': 12},
 '2ACY': {'A': 98},
 '2ADA': {'A': 349},
 '2AFG': {'A': 129, 'B': 129, 'C': 127, 'D': 128},
 '2AKY': {'A': 218},
 '2BQA': {'A': 130},
 '2BRD': {'A': 222},
 '2CI2': {'I': 65},
 '2CRK': {'A': 365},
 '2HMB': {'A': 131},
 '2HPR': {'A': 87},
 '2IFB': {'A': 131},
 '2IMM': {'A': 114},
 '2L3Y': {'A': 184},
 '2LZM': {'A': 164},
 '2Q98': {'A': 191},
 '2RN2': {'A': 155},
 '2TRX': {'A': 108, 'B': 108},
 '2TS1': {'A': 317},
 '2WSY': {'A': 237, 'B': 375},
 '2ZTA': {'A': 31, 'B': 31},
 '3BCI': {'A': 165},
 '3BLS': {'A': 357, 'B': 357},
 '3D2A': {'A': 179},
 '3HHR': {'A': 185, 'B': 197, 'C': 196},
 '3K0NA_lin': {'A': 163},
 '3K0NB_lin': {'A': 163},
 '3K0On_lin': {'A': 164},
 '3MBP': {'A': 370},
 '3PGK': {'A': 415},
 '3SSI': {'A': 108},
 '451C': {'A': 82},
 '4BLM': {'A': 256, 'B': 256},
 '4LYZ': {'A': 129},
 '5AZU': {'A': 128, 'B': 128, 'C': 128, 'D': 128},
 '5CRO': {'A': 61, 'B': 61, 'C': 61, 'O': 60},
 'S9G10_best': {'A': 157, 'B': 369},
 'ub_CUE': {'A': 76, 'B': 52},
 'ub_OTU': {'A': 75, 'B': 172},
 'ub_RPN13': {'A': 76, 'B': 109},
 'ub_SH3': {'A': 76, 'B': 71},
 'ub_UQcon': {'A': 76, 'B': 152, 'C': 135},
 'uby_1UBQ': {'A': 76},
 'uby_CUE': {'A': 76, 'B': 52},
 'uby_OTU': {'A': 75, 'B': 172},
 'uby_RPN13': {'A': 76, 'B': 109},
 'uby_SH3': {'A': 76, 'B': 71},
 'uby_UQcon': {'A': 76, 'B': 152, 'C': 135}}

def analyze_data():
    ddGdb = ddgdbapi.ddGDatabase()

    global pdb_chain_lengths
    if not pdb_chain_lengths:
        pdb_chain_lengths = {}
        pdb_os = {}
        for r in ddGdb.execute_select('SELECT DISTINCT PDBFileID FROM Prediction INNER JOIN Experiment ON ExperimentID=Experiment.ID WHERE Status="done" AND maxvmem IS NOT NULL'):
            pdb_id = r['PDBFileID']
            pdb_chain_lengths[pdb_id] = pdb_chain_lengths.get(pdb_id, {})
            p = PDB(ddGdb.execute_select('SELECT Content FROM PDBFile WHERE ID=%s', parameters=(pdb_id))[0]['Content'])
            for chain_id, s in p.atom_sequences.iteritems():
                pdb_chain_lengths[pdb_id][chain_id] = len(s)

        import pprint
        pprint.pprint(pdb_chain_lengths)


    # Gather data points
    mem_points = []
    time_points = []
    for r in ddGdb.execute_select('SELECT Prediction.ID, ExperimentID, maxvmem, DDGTime, PDBFileID FROM Prediction INNER JOIN Experiment ON Experiment.ID=ExperimentID WHERE Status="done" AND PredictionSet="Protocol_16_r57471" ORDER BY ID'):
        PredictionID = r['ID']
        ExperimentID = r['ExperimentID']
        maxvmem = r['maxvmem']
        DDGTime = r['DDGTime']
        pdb_id = r['PDBFileID']
        chains = [cr for cr in ddGdb.execute_select('SELECT * FROM ExperimentChain WHERE ExperimentID=%s', parameters=(ExperimentID,))]
        assert(len(chains) == 1)
        chain = chains[0]['Chain']
        if maxvmem != None:
            mem_points.append((pdb_chain_lengths[pdb_id][chain], maxvmem, PredictionID))
        if DDGTime != None:
            time_points.append((pdb_chain_lengths[pdb_id][chain], DDGTime, PredictionID))

    # Create input files
    lines = ['ChainLength,MaxVMem,PredictionID']
    for p in mem_points:
        lines.append('%s,%s,%s' % (str(p[0]), str(p[1]), str(p[2])))
    write_file('extract_maxvmem.txt', '\n'.join(lines))
    lines = ['ChainLength,DDGTime,PredictionID']
    for p in time_points:
        lines.append('%s,%s,%s' % (str(p[0]), str(p[1]), str(p[2])))
        print(lines[-1])
    write_file('extract_ddgtime.txt', '\n'.join(lines))

    # Determine limits
    success = False
    min_constant = 1.2 # set a minimum memory cap
    max_scalar = float(300) # the amount of requested memory is divided by this number
    failed_cases = 0
    while max_scalar > 0:
        #colortext.message('Trying max_scalar=%0.2f.' % max_scalar)
        failed_cases = 0
        for p in mem_points:
            est_max = max(min_constant, float(p[0]) / max_scalar)
            actual_max = p[1]
            if est_max < actual_max:
                failed_cases += 1
        if failed_cases <= 10:
            success = True
            break
        max_scalar -= 10

    if success:
        print('Setting the memory requirement to max(%0.2f, len(chain)/%d) GB  would have been sufficient '
              'for all but %d of the %d runs in this prediction set.' % (min_constant, max_scalar, failed_cases, len(mem_points)))
        write_file("extract_maxvmem_temp.R", read_file("extract_maxvmem.R").replace('LIMITNOTE', 'Jobs require max(%0.2f, len(chain)/%d) GB' % (min_constant, max_scalar)))

        # Create an output file for R
        p = subprocess.Popen(["R", "CMD", "BATCH", "extract_maxvmem_temp.R"])
        while True:
            time.sleep(0.3)
            errcode = p.poll()
            if errcode != None:
                break
        os.remove("extract_maxvmem_temp.R")
        os.remove("extract_maxvmem_temp.Rout")


    # Determine limits
    success = False
    min_constant = 300.0 # set a minimum time cap
    max_scalar = float(0.2) # the amount of requested time is divided by this number
    failed_cases = 0
    while max_scalar > 0:
        failed_cases = 0
        for p in time_points:
            est_max = max(min_constant, float(p[0]) / max_scalar)
            actual_max = p[1]
            if est_max < actual_max:
                failed_cases += 1
        if failed_cases == 0: # hard limit
            success = True
            break
        max_scalar -= 0.0001

    if success:
        for p in time_points:
            print(p)
            if p[0] <= 300:
                assert(2000 > p[1]) # 2400 = 40 hours should be safe => set to 2 days
            elif p[0] <= 400:
                assert(3000 > p[1]) # similarly, set to 3 days
            elif p[0] <= 500:
                assert(5000 > p[1]) # this is 3.5 days so set the max time to 1 week
            else:
                assert(20000 > p[1]) # set the max time to 2 weeks
            est_max = max(min_constant, float(p[0]) / max_scalar)
            actual_max = p[1]

    if success:
        print('Setting the time requirement to max(%0.2f, len(chain)/%f) hours would have been sufficient '
              'for all %d successful runs in this prediction set.' % (min_constant/60.0, max_scalar*60, len(time_points)))
        write_file("extract_maxtime_temp.R", read_file("extract_maxtime.R").replace('LIMITNOTE', 'Jobs require max(%0.2f, len(chain)/%f) hours' % (min_constant/60.0, max_scalar*60)))

        # Create an output file for R
        p = subprocess.Popen(["R", "CMD", "BATCH", "extract_maxtime_temp.R"])
        while True:
            time.sleep(0.3)
            errcode = p.poll()
            if errcode != None:
                break
        os.remove("extract_maxtime_temp.R")
        os.remove("extract_maxtime_temp.Rout")
        print('Done.')


if __name__ == '__main__':
    #extract_data()
    analyze_data()