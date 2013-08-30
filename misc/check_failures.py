import sys
import os

sys.path.insert(0, "../..")
sys.path.insert(0, "..")
from tools import colortext
from tools.fs.io import read_file
from ddglib import ddgdbapi

def ask_yes_no(question_text, default_value = None, max_tries = None):
    '''Returns true if the user answers 'y' or 'yes' and False if the user answers 'n' or 'no'. If the user does not
    provide an answer and default_value is not None, default_value is returned as the answer. If the answer is invalid,
    the user will be repeatedly asked. If max_tries is set and exceeded then None is returned. '''
    answer = None
    num_tries = 0
    while answer == None:
        if question_text:
            print("%s " % question_text)
        a = raw_input().strip().lower()
        if a == 'y' or a == 'yes':
            answer = True
        elif a == 'n' or a == 'no':
            answer = False
        elif a == '' and default_value != None:
            answer = default_value
        num_tries += 1
        if max_tries and num_tries >= max_tries:
            return None
    return answer

def check_failures(prediction_set):
    ddGdb = ddgdbapi.ddGDatabase()
    results_root = '/var/cluster/temp/'

    UserDataSetExperimentIDs = {}
    results = ddGdb.execute_select('''SELECT ID, ExperimentID FROM Prediction WHERE PredictionSet = %s AND STATUS = 'failed' ''', parameters=(prediction_set,))
    reported_failures = [r['ID'] for r in results]
    for r in results:
        UserDataSetExperimentIDs[r['ID']] = r['ExperimentID']

    affected_subsets = {}
    actually_failed = []
    did_not_fail = []
    for PredictionID in reported_failures:
        print(PredictionID)
        results_dir = os.path.join(results_root, prediction_set, str(PredictionID))
        print(results_dir)
        file_list = os.listdir(results_dir)
        print(file_list)

        found_stdout = 0
        found_stderr = 0
        for f in file_list:
            if f.find('.cmd.o') != -1:
                found_stdout = 1
            elif f.find('.cmd.e') != -1:
                found_stderr = 1
                colortext.error(f)
        assert(found_stdout >= found_stderr)
        if found_stderr:
            assert(found_stderr == 1)
            colortext.error("Job #%d actually failed" % PredictionID)
            actually_failed.append(PredictionID)

            ExperimentID = UserDataSetExperimentIDs[PredictionID]
            sub_results = ddGdb.execute_select('''SELECT Subset FROM UserAnalysisSet WHERE ExperimentID = %s''', parameters=(ExperimentID,))
            for sr in sub_results:
                print(sr['Subset'])

                affected_subsets[sr['Subset']] = affected_subsets.get(sr['Subset'], [])
                affected_subsets[sr['Subset']].append(PredictionID)



        else:
            colortext.warning("Job #%d had not failed by the time it was terminated." % PredictionID)
            did_not_fail.append(PredictionID)

    colortext.message("*** Report ***")
    print('%d jobs were marked as failed.' % len(reported_failures))
    colortext.warning('%d jobs were marked as failed but had not failed.' % len(did_not_fail))
    colortext.error('%d jobs were marked as failed and did fail.\n' % len(actually_failed))

    if affected_subsets:
        print("The following subsets were affected")
        for k, v in affected_subsets.iteritems():
            print("%s: %d records" % (k, len(v)))

    restart_jobs = ask_yes_no("Do you want to restart the jobs that did not actually fail?", default_value=False)
    if restart_jobs:
        for j in did_not_fail:
            r = ddGdb.execute("SELECT Status, AdminCommand FROM Prediction WHERE ID=%s", parameters=(j,))
            assert(len(r) == 1)
            if r[0]['Status'] != 'failed' or r[0]['AdminCommand'] != 'restart':
                ddGdb.execute("UPDATE Prediction SET AdminCommand='restart' WHERE ID=%s", parameters=(j,))

    for k, v in affected_subsets.iteritems():
        see_errors = ask_yes_no("Do you want to see the stderr files for the jobs that did fail and affected %s?" % k, default_value=False)
        if see_errors:
            count = 1
            for j in v:
                results_dir = os.path.join(results_root, prediction_set, str(j))
                file_list = os.listdir(results_dir)

                colortext.error("\n[%d/%d] Prediction ID: %d" % (count, len(v), j))
                PDBFileID = ddGdb.execute_select("SELECT UserDataSetExperiment.PDBFileID AS PDBFileID FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperimentID=UserDataSetExperiment.ID WHERE Prediction.ID=%s", parameters=(j,))
                assert(len(PDBFileID) == 1)
                PDBFileID = PDBFileID[0]['PDBFileID']
                colortext.error("\nPDB ID: %s" % PDBFileID)
                for f in file_list:
                    if f.find('.cmd.e') != -1:
                        colortext.warning(f)
                        print(read_file(os.path.join(results_dir, f)))
                        print("")
                count += 1

def fix_1TEN_InputFiles():
    '''This is a once-off function which should only be run once as otherwise'''
    import pickle
    ddGdb = ddgdbapi.ddGDatabase()

    BadPredictions = sorted(set([(r['PredictionID'], r['Status']) for r in ddGdb.execute_select("SELECT Prediction.ID AS PredictionID, Status FROM Prediction INNER JOIN Experiment ON Experiment.ID=Prediction.ExperimentID INNER JOIN ExperimentMutation ON Experiment.ID=ExperimentMutation.ExperimentID WHERE PredictionSet = 'RosCon2013_P16_talaris2013' AND PDBFileID='1TEN' ")]))
    BadPredictionIDs = sorted(set([r[0] for r in BadPredictions]))
    num_active = len([r for r in BadPredictions if r[1] == 'active'])
    num_queued = len([r for r in BadPredictions if r[1] == 'queued'])
    statuses = sorted(set([r[1] for r in BadPredictions]))
    if ('active' in statuses) or ('queued' in statuses):
        colortext.error("Cannot proceed - there are %d active jobs and %d queued in the list that need to be fixed up. Stop the DDG scheduler, remove the queued constraint, and rerun this function. " % (num_active, num_queued))
        if num_active:
            print("%d active jobs: %s" % (num_active, ", ".join([str(r[0]) for r in BadPredictions if r[1] == 'active'])))
        if num_queued:
            print("%d queued jobs: %s" % (num_queued, ", ".join([str(r[0]) for r in BadPredictions if r[1] == 'queued'])))
        #3return

    for PredictionID in BadPredictionIDs:
        r = ddGdb.execute_select("SELECT InputFiles FROM Prediction WHERE ID=%s", parameters=(PredictionID,))
        assert(len(r) == 1)
        r = r[0]

        InputFiles = pickle.loads(r['InputFiles'])
        assert(InputFiles.keys() == ['MUTFILE'])
        mutfile = InputFiles['MUTFILE']

        colortext.message("\n%d" % PredictionID)

        colortext.warning('original')
        print(mutfile)

        lines = mutfile.split("\n")
        assert(lines[0].startswith('total'))
        num_muts = int(lines[0][5:])
        assert(lines[1] == str(num_muts))
        for x in range(2, num_muts + 2):
            mutline = lines[x]
            tokens = mutline.split()
            tokens[1] = str(int(tokens[1]) - 1)
            lines[x] = " ".join(tokens)

        new_mutfile = "\n".join(lines)
        colortext.warning('fixed')
        print(new_mutfile)

        p = pickle.dumps({'MUTFILE' : new_mutfile})

        #results = ddGdb.execute("UPDATE Prediction SET Status='queued', InputFiles=%s WHERE ID=%s", parameters=(p, PredictionID))

def count_num_residues_in_active_jobs():
    '''I wrote this function to try to narrow down into which jobs ran the longest as I suspect that this is due to long PDB chains.'''
    ddGdb = ddgdbapi.ddGDatabase()
    active_jobs = ddGdb.execute_select("SELECT DISTINCT ExperimentID FROM Prediction WHERE Status='active'")
    colortext.message("\n%d jobs are active" % len(active_jobs))

    from tools.bio.rcsb import parseFASTAs

    chains_in_active_jobs = {}
    PDB_chain_lengths ={}
    for active_job in active_jobs:
        r = ddGdb.execute_select('SELECT PDBFileID, Chain FROM Experiment INNER JOIN ExperimentChain ON ExperimentID=Experiment.ID WHERE ExperimentID=%s', parameters=(active_job['ExperimentID']))
        assert(len(r) == 1)
        r = r[0]

        key = (r['PDBFileID'], r['Chain'])

        if PDB_chain_lengths.get(key) == None:
            fasta = ddGdb.execute_select("SELECT FASTA FROM PDBFile WHERE ID=%s", parameters = (r['PDBFileID'],))
            assert(len(fasta) == 1)
            fasta = fasta[0]['FASTA']
            f = parseFASTAs(fasta)
            PDB_chain_lengths[key] = len(f[r['PDBFileID']][r['Chain']])
        chain_length = PDB_chain_lengths[key]

        chains_in_active_jobs[key] = chains_in_active_jobs.get(key, [chain_length, 0])
        chains_in_active_jobs[key][1] += 1

    if chains_in_active_jobs:
        colortext.message("Chains in currently active jobs:\n")

        print("PDB\tChain\tChain SEQRES length\tJobs remaining")
        for k,v in sorted(chains_in_active_jobs.iteritems(), key=lambda x: x[1][0]):
            print("%s\t  %s\t%s\t%s" % (k[0], k[1], str(v[0]).center(19), str(v[1]).center(14)))


#check_failures('RosettaCon2013_P16_talaris2013')

check_failures('RosCon2013_P16_score12prime')
#check_failures('RosCon2013_P16_talaris2013')
#fix_1TEN_InputFiles()
#count_num_residues_in_active_jobs()
