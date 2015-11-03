import sys
import os
import zipfile

sys.path.insert(0, "../..")
sys.path.insert(0, "..")
from tools import colortext
from tools.fs.fsio import read_file
from tools.bio.pdb import PDB
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
    results_root = '/kortemmelab/shared/DDG/jobs'

    UserDataSetExperimentIDs = {}
    results = ddGdb.execute_select('''SELECT ID, ExperimentID FROM Prediction WHERE PredictionSet = %s AND STATUS = 'failed' ''', parameters=(prediction_set,))
    reported_failures = [r['ID'] for r in results]
    for r in results:
        UserDataSetExperimentIDs[r['ID']] = r['ExperimentID']
    print(UserDataSetExperimentIDs)
    print(len(UserDataSetExperimentIDs))

    affected_subsets = {}
    actually_failed = []
    did_not_fail = []
    did_not_fail_but_have_many_residues = []
    did_not_fail_but_has_a_troublesome_structure = []
    large_proteins = set(['1FEP', '1W99'])
    troublesome_structures = set(['2IFB', '1FMK'])

    for PredictionID in reported_failures:
        print(PredictionID)
        zipfile_path = os.path.join(results_root, '%d.zip' % PredictionID)
        #try:
        z = zipfile.ZipFile(zipfile_path, 'r')
        #except:
        #    colortext.error('MISSING FILE FOR %d' % PredictionID)
        #    did_not_fail.append(PredictionID)
        #    continue
        file_list = z.namelist()

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
            PDBFileID = ddGdb.execute_select("SELECT UserDataSetExperiment.PDBFileID AS PDBFileID FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperimentID=UserDataSetExperiment.ID WHERE Prediction.ID=%s", parameters=(PredictionID,))[0]['PDBFileID']
            if PDBFileID in troublesome_structures:
                colortext.warning("Job #%d had not failed by the time it was terminated however it has a troublesome structure (%s)." % (PredictionID, PDBFileID))
                did_not_fail_but_has_a_troublesome_structure.append(PredictionID)
            elif PDBFileID in large_proteins:
                colortext.warning("Job #%d had not failed by the time it was terminated however it has many residues (%s)." % (PredictionID, PDBFileID))
                did_not_fail_but_have_many_residues.append(PredictionID)
            else:
                colortext.warning("Job #%d had not failed by the time it was terminated." % PredictionID)
                did_not_fail.append(PredictionID)


    colortext.message("*** Report ***")
    print('%d jobs were marked as failed.' % len(reported_failures))
    colortext.warning('%d jobs were marked as failed but had not failed.' % len(did_not_fail))
    colortext.warning('%d jobs were marked as failed and had not failed but had large chains.' % len(did_not_fail_but_have_many_residues))
    colortext.warning('%d jobs were marked as failed and had not failed but have a troublesome structure.' % len(did_not_fail_but_has_a_troublesome_structure))
    colortext.error('%d jobs were marked as failed and did fail.\n' % len(actually_failed))

    if affected_subsets:
        print("The following subsets were affected")
        for k, v in affected_subsets.iteritems():
            print("%s: %d records" % (k, len(v)))

    if did_not_fail:
        restart_jobs = ask_yes_no("Do you want to restart the jobs that did not actually fail?", default_value=False)
        if restart_jobs:
            for j in did_not_fail:
                r = ddGdb.execute("SELECT Status, AdminCommand FROM Prediction WHERE ID=%s", parameters=(j,))
                assert(len(r) == 1)
                if r[0]['Status'] != 'failed' or r[0]['AdminCommand'] != 'restart':
                    ddGdb.execute("UPDATE Prediction SET AdminCommand='restart' WHERE ID=%s", parameters=(j,))

    if did_not_fail_but_have_many_residues:
        restart_jobs = ask_yes_no("Do you want to restart the jobs that did not actually fail but had large proteins?", default_value=False)
        if restart_jobs:
            for j in did_not_fail_but_have_many_residues:
                r = ddGdb.execute("SELECT Status, AdminCommand FROM Prediction WHERE ID=%s", parameters=(j,))
                assert(len(r) == 1)
                if r[0]['Status'] != 'failed' or r[0]['AdminCommand'] != 'restart':
                    ddGdb.execute("UPDATE Prediction SET AdminCommand='restart' WHERE ID=%s", parameters=(j,))

    if did_not_fail_but_has_a_troublesome_structure:
        restart_jobs = ask_yes_no("Do you want to restart the jobs that did not actually fail but had troublesome structures?", default_value=False)
        if restart_jobs:
            for j in did_not_fail_but_has_a_troublesome_structure:
                r = ddGdb.execute("SELECT Status, AdminCommand FROM Prediction WHERE ID=%s", parameters=(j,))
                assert(len(r) == 1)
                if r[0]['Status'] != 'failed' or r[0]['AdminCommand'] != 'restart':
                    ddGdb.execute("UPDATE Prediction SET AdminCommand='restart' WHERE ID=%s", parameters=(j,))

    for k, v in affected_subsets.iteritems():
        see_errors = ask_yes_no("Do you want to see the stderr files for the jobs that did fail and affected %s?" % k, default_value=False)
        if see_errors:
            count = 1
            for j in v:
                zipfile_path = os.path.join(results_root, '%d.zip' % j)
                z = zipfile.ZipFile(zipfile_path, 'r')
                file_list = z.namelist()

                colortext.error("\n[%d/%d] Prediction ID: %d" % (count, len(v), j))
                PDBFileID = ddGdb.execute_select("SELECT UserDataSetExperiment.PDBFileID AS PDBFileID FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperimentID=UserDataSetExperiment.ID WHERE Prediction.ID=%s", parameters=(j,))
                assert(len(PDBFileID) == 1)
                PDBFileID = PDBFileID[0]['PDBFileID']
                colortext.error("\nPDB ID: %s" % PDBFileID)
                for f in file_list:
                    if f.find('.cmd.e') != -1:
                        colortext.warning(f)
                        print(z.open(f, 'r').read()[:300])
                        print("")
                count += 1

def fix_1TEN_InputFiles(prediction_set):
    '''This is a once-off function which should only be run once per prediction set as each run changes the mutfile and this change should only occur once.'''
    import pickle
    ddGdb = ddgdbapi.ddGDatabase()

    BadPredictions = sorted(set([(r['PredictionID'], r['Status']) for r in ddGdb.execute_select("SELECT Prediction.ID AS PredictionID, Status FROM Prediction INNER JOIN Experiment ON Experiment.ID=Prediction.ExperimentID INNER JOIN ExperimentMutation ON Experiment.ID=ExperimentMutation.ExperimentID WHERE PredictionSet=%s AND PDBFileID='1TEN' ", parameters=(prediction_set,))]))
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
        return

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

def fix_1H7M_InputFiles(prediction_set):
    '''This is a once-off function which should only be run once per prediction set as each run changes the mutfile and this change should only occur once.'''
    import pickle
    ddGdb = ddgdbapi.ddGDatabase()

    BadPredictions = sorted(set([(r['PredictionID'], r['Status']) for r in ddGdb.execute_select('''
    SELECT Prediction.ID AS PredictionID, Status FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperiment.ID=Prediction.UserDataSetExperimentID WHERE PredictionSet=%s AND PDBFileID='1H7M'
    ''', parameters=(prediction_set,))]))
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
        return

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

def fix_1AYE_InputFiles(prediction_set):
    '''This is a once-off function which should only be run once per prediction set as each run changes the mutfile and this change should only occur once.'''
    import pickle
    ddGdb = ddgdbapi.ddGDatabase()

    BadPredictions = sorted(set([(r['PredictionID'], r['Status']) for r in ddGdb.execute_select('''
    SELECT Prediction.ID AS PredictionID, Status FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperiment.ID=Prediction.UserDataSetExperimentID WHERE PredictionSet=%s AND PDBFileID='1AYE'
    ''', parameters=(prediction_set,))]))
    BadPredictionIDs = sorted(set([r[0] for r in BadPredictions]))
    print(BadPredictions)
    num_active = len([r for r in BadPredictions if r[1] == 'active'])
    num_queued = len([r for r in BadPredictions if r[1] == 'queued'])
    statuses = sorted(set([r[1] for r in BadPredictions]))
    if ('active' in statuses) or ('queued' in statuses):
        colortext.error("Cannot proceed - there are %d active jobs and %d queued in the list that need to be fixed up. Stop the DDG scheduler, remove the queued constraint, and rerun this function. " % (num_active, num_queued))
        if num_active:
            print("%d active jobs: %s" % (num_active, ", ".join([str(r[0]) for r in BadPredictions if r[1] == 'active'])))
        if num_queued:
            print("%d queued jobs: %s" % (num_queued, ", ".join([str(r[0]) for r in BadPredictions if r[1] == 'queued'])))
        return

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

def remove_existing_zip_files(prediction_set):
    raise Exception('This is a clean-up function. You probably do not want to run for fixing runs this since it interferes with check_failures above.')
    ddGdb = ddgdbapi.ddGDatabase()
    results_root = '/kortemmelab/shared/DDG/jobs'
    ids = ddGdb.execute_select("SELECT ID FROM Prediction WHERE Status IN ('failed', 'queued', 'active') AND PredictionSet=%s", parameters=(prediction_set,))
    #ids = ddGdb.execute_select("SELECT ID FROM Prediction WHERE Status IN ('done') AND PredictionSet=%s", parameters=(prediction_set,))
    count = 0
    for id in ids:
        if os.path.exists(os.path.join(results_root, '%s.zip' % id['ID'])):
            os.remove(os.path.join(results_root, '%s.zip' % id['ID']))
            count += 1
    print(len(ids))
    print(count)

def classify_failures(prediction_set):
    ddGdb = ddgdbapi.ddGDatabase()
    results_root = '/kortemmelab/shared/DDG/jobs'

    UserDataSetExperimentIDs = {}
    results = ddGdb.execute_select('''SELECT ID, ExperimentID FROM Prediction WHERE PredictionSet = %s AND STATUS = 'failed' ''', parameters=(prediction_set,))
    reported_failures = [r['ID'] for r in results]
    for r in results:
        UserDataSetExperimentIDs[r['ID']] = r['ExperimentID']

    actually_failed = []
    did_not_fail = []
    for PredictionID in reported_failures:
        zipfile_path = os.path.join(results_root, '%d.zip' % PredictionID)
        #try:
        z = zipfile.ZipFile(zipfile_path, 'r')
        #except:
        #    colortext.error('MISSING FILE FOR %d' % PredictionID)
        #    continue
        file_list = z.namelist()

        found_stdout = 0
        found_stderr = 0
        for f in file_list:
            if f.find('.cmd.o') != -1:
                found_stdout = 1
            elif f.find('.cmd.e') != -1:
                found_stderr = 1
        assert(found_stdout >= found_stderr)

        if found_stderr:
            assert(found_stderr == 1)
            colortext.error("Job #%d actually failed" % PredictionID)
            actually_failed.append(PredictionID)
        else:
            colortext.warning("Job #%d had not failed by the time it was terminated." % PredictionID)
            did_not_fail.append(PredictionID)

    colortext.message("*** Report ***")
    print('%d jobs were marked as failed.' % len(reported_failures))
    colortext.warning('%d jobs were marked as failed but had not failed.' % len(did_not_fail))
    colortext.error('%d jobs were marked as failed and did fail.\n' % len(actually_failed))

    pdb_details = {}
    failed_job_pdb_files = {}
    for failed_job in actually_failed:
        PDBFileID = ddGdb.execute_select("SELECT UserDataSetExperiment.PDBFileID AS PDBFileID FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperimentID=UserDataSetExperiment.ID WHERE Prediction.ID=%s", parameters=(failed_job,))[0]['PDBFileID']
        pdb_details[PDBFileID] = True
        failed_job_pdb_files[failed_job] = PDBFileID

    for pdb_id in pdb_details.keys():
        pdb_details[pdb_id] = ddGdb.execute_select("SELECT Resolution, Techniques FROM PDBFile WHERE ID=%s", parameters=(pdb_id,))[0]
        pdb_details[pdb_id]['Chains'] = [r['Chain'] for r in ddGdb.execute_select("SELECT Chain FROM PDBChain WHERE PDBFileID=%s ORDER BY Chain", parameters=(pdb_id,))]
        pdb_details[pdb_id]['TotalJobs'] = ddGdb.execute_select("SELECT Count(ID) AS TotalJobs FROM UserDataSetExperiment WHERE PDBFileID=%s", parameters=(pdb_id,))[0]['TotalJobs']

    hosts = {}
    failed_by_hessin = {}
    failed_by_residue_mismatch = {}
    failed_for_another_reason = {}
    missing_output = {}
    mutfiles = {}
    count = 1
    for failed_job in actually_failed:
        mutfile = None
        colortext.message('Failed job %d of %d' % (count, failed_job))
        zipfile_path = os.path.join(results_root, '%d.zip' % failed_job)
        found_output = False
        pdb_id = failed_job_pdb_files[failed_job]
        if os.path.exists(zipfile_path):
            z = zipfile.ZipFile(zipfile_path, 'r')
            file_list = z.namelist()
            for f in file_list:
                if f.find('.cmd.e') != -1:
                    found_output = True
                    stderr_contents = z.open(f, 'r').read()
                    stdout_contents = z.open(f.replace('.cmd.e', '.cmd.o'), 'r').read()

                    hosts[failed_job] = stdout_contents[stdout_contents.find('<host>') + 6:stdout_contents.find('</host>')].strip()

                    if stderr_contents.find('HESSIN for (i,i):') != -1:
                        assert(stderr_contents.find('G for (i):') != -1)
                        print(stderr_contents[:120])
                        failed_by_hessin[pdb_id] = failed_by_hessin.get(pdb_id, [])
                        failed_by_hessin[pdb_id].append(failed_job)
                        colortext.error('HESSIN: %s' % pdb_id)
                    elif stderr_contents.find('ERROR: pose.residue(resnum).name1() == wt') != -1:
                        failed_by_residue_mismatch[pdb_id] = failed_by_residue_mismatch.get(pdb_id, [])
                        failed_by_residue_mismatch[pdb_id].append(failed_job)
                        colortext.error('MISMATCH')
                    else:
                        failed_for_another_reason[pdb_id] = failed_for_another_reason.get(pdb_id, [])
                        failed_for_another_reason[pdb_id].append(failed_job)
                        colortext.error('UNKNOWN')
                        see_errors = ask_yes_no("Do you want to see the stderr files for prediction %d?" % failed_job, default_value=False)
                        if see_errors:
                            colortext.warning(f)
                            print(stderr_contents[:300])
                            print("")
                if f.find('.mutfile') != -1:
                    assert(mutfile == None)
                    mutfile = z.open(f, 'r').read()
                    mutfiles[failed_job] = mutfile

        if not found_output:
            missing_output[pdb_id] = missing_output.get(pdb_id, [])
            missing_output[pdb_id].append(failed_job)
        count += 1


    colortext.message("*** Report ***")
    if missing_output:
        colortext.warning("Missing output: %d jobs" % sum([len(v) for k, v in missing_output.iteritems()]))
        for k, v in sorted(missing_output.iteritems()):
            print('%s: %d jobs - %s' % (k, len(v), ', '.join(map(str, sorted(v)))))
    if failed_by_hessin:
        colortext.warning("Failed Hessin: %d jobs" % sum([len(v) for k, v in failed_by_hessin.iteritems()]))
        for k, v in sorted(failed_by_hessin.iteritems()):
            if pdb_details[k]['Resolution'] != None:
                print('%s, %0.2fA, %s.' % (k, pdb_details[k]['Resolution'], pdb_details[k]['Techniques'].title()))
            else:
                print('%s, %s.' % (k, pdb_details[k]['Techniques'].title()))
            print('%d/%d jobs failed - %s\n' % (len(v), pdb_details[k]['TotalJobs'], ', '.join(map(str, sorted(v)))))
            for failed_id in sorted(v):
                mutations = ddGdb.execute_select("SELECT Prediction.ExperimentID, ExperimentMutation.* FROM Prediction INNER JOIN ExperimentMutation ON Prediction.ExperimentID=ExperimentMutation.ExperimentID WHERE Prediction.ID=%s", parameters=(failed_id,))
                mut_str = ', '.join([('%s %s%s%s' % (m['Chain'], m['WildTypeAA'], m['ResidueID'], m['MutantAA'])) for m in mutations])
                colortext.printf('%d: %s, experiment #%d. Host = %s' % (failed_id, mut_str, mutations[0]['ExperimentID'], hosts[failed_id]), 'orange')
                print('')
            print('')

    if failed_by_residue_mismatch:
        colortext.warning("Failed due to residue mismatch: %d jobs" % sum([len(v) for k, v in failed_by_residue_mismatch.iteritems()]))
        for k, v in sorted(failed_by_residue_mismatch.iteritems()):
            if pdb_details[k]['Resolution'] != None:
                colortext.printf('%s, %0.2fA, %s.' % (k, pdb_details[k]['Resolution'], pdb_details[k]['Techniques'].title()), 'cyan')
            else:
                colortext.printf('%s, %s.' % (k, pdb_details[k]['Techniques'].title()), 'cyan')
            print('%d/%d jobs failed - %s\n' % (len(v), pdb_details[k]['TotalJobs'], ', '.join(map(str, sorted(v)))))
            for failed_id in sorted(v):
                mutations = ddGdb.execute_select("SELECT ExperimentMutation.* FROM Prediction INNER JOIN ExperimentMutation ON Prediction.ExperimentID=ExperimentMutation.ExperimentID WHERE Prediction.ID=%s", parameters=(failed_id,))
                mut_str = ', '.join([('%s %s%s%s' % (m['Chain'], m['WildTypeAA'], m['ResidueID'], m['MutantAA'])) for m in mutations])
                colortext.printf('%d: %s' % (failed_id, mut_str), 'orange')
                print(mutfiles[failed_id])
                print('')
            print('')

    if failed_for_another_reason:
        colortext.warning("Failed for an unknown reason: %d jobs" % sum([len(v) for k, v in failed_for_another_reason.iteritems()]))
        for k, v in sorted(failed_for_another_reason.iteritems()):
            if pdb_details[k]['Resolution'] != None:
                print('%s, %0.2fA, %s.' % (k, pdb_details[k]['Resolution'], pdb_details[k]['Techniques'].title()))
            else:
                print('%s, %s.' % (k, pdb_details[k]['Techniques'].title()))
            print('%d/%d jobs failed - %s\n' % (len(v), pdb_details[k]['TotalJobs'], ', '.join(map(str, sorted(v)))))

    print('%d jobs were marked as failed.' % len(reported_failures))
    colortext.warning('%d jobs were marked as failed but had not failed.' % len(did_not_fail))
    colortext.error('%d jobs were marked as failed and did fail.\n' % len(actually_failed))


pdb_chain_lengths = {
 '1A23': {'A': 189},
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
 'uby_UQcon': {'A': 76, 'B': 152, 'C': 135}
}

def add_job_costs(prediction_set):
    ddGdb = ddgdbapi.ddGDatabase()
    global pdb_chain_lengths
    if not pdb_chain_lengths:
        pdb_chain_lengths = {}
        pdb_os = {}
        for r in ddGdb.execute_select('SELECT DISTINCT PDBFileID FROM Prediction INNER JOIN Experiment ON ExperimentID=Experiment.ID'):
            pdb_id = r['PDBFileID']
            pdb_chain_lengths[pdb_id] = pdb_chain_lengths.get(pdb_id, {})
            p = PDB(ddGdb.execute_select('SELECT Content FROM PDBFile WHERE ID=%s', parameters=(pdb_id))[0]['Content'])
            for chain_id, s in p.atom_sequences.iteritems():
                pdb_chain_lengths[pdb_id][chain_id] = len(s)

        import pprint
        pprint.pprint(pdb_chain_lengths)

    for r in ddGdb.execute_select('SELECT Prediction.ID, Status, PDBFileID, Cost FROM Prediction INNER JOIN Experiment ON Experiment.ID=ExperimentID WHERE PredictionSet=%s AND Status<>"done" ORDER BY Cost', parameters=(prediction_set,)):
        print(r)
    return
    for r in ddGdb.execute_select('SELECT Prediction.ID, ExperimentID, PDBFileID FROM Prediction INNER JOIN Experiment ON Experiment.ID=ExperimentID WHERE PredictionSet=%s ORDER BY ID', parameters=(prediction_set,)):
        PredictionID = r['ID']
        ExperimentID = r['ExperimentID']
        pdb_id = r['PDBFileID']
        chains = [cr for cr in ddGdb.execute_select('SELECT * FROM ExperimentChain WHERE ExperimentID=%s', parameters=(ExperimentID,))]
        assert(len(chains) == 1)
        chain = chains[0]['Chain']
        cost = pdb_chain_lengths[pdb_id][chain]
        ddGdb.execute('UPDATE Prediction SET Cost=%s WHERE ID=%s', parameters=(cost, PredictionID,))


#remove_existing_zip_files('Protocol_16_r57471')

check_failures('Protocol_16_r57471')
#add_job_costs('Protocol_16_r57471')
#classify_failures('Protocol_16_r57471')

#fix_1H7M_InputFiles('Protocol_16_r57471')
#fix_1AYE_InputFiles('Protocol_16_r57471')
#fix_1TEN_InputFiles('Protocol_16_r57471')

#check_failures('RosCon2013_P16_score12prime')
#check_failures('RosCon2013_P16_talaris2013')
#fix_1TEN_InputFiles()
#count_num_residues_in_active_jobs()
