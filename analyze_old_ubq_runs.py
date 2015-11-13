import os, sys
import shutil
import pprint

import klab.cluster_template.parse_settings as parse_settings
from klab.Reporter import Reporter
import time
import getpass
import tempfile
import json
import re
import datetime
import subprocess
from klab.cluster_template.write_run_file import process as write_run_file
from ddglib.ddg_monomer_ppi_api import get_interface as get_interface_factory
from ddglib.ppi_api import get_interface_with_config_file
from klab.fs.zip_util import unzip_file

from klab import colortext
from klab.bio.pdb import PDB

from ddglib import ddgdbapi, db_api

import cPickle as pickle

tmpdir_location = '/dbscratch/%s/tmp' % getpass.getuser()

ubiquitin_chains = [
    ('1ubq', 'A', 'Ubiquitin scan: 1UBQ p16'),
    ('uby_1UBQ', 'A', 'Ubiquitin scan: 1UBQ_yeast p16'),
    ('ub_SH3', 'A', 'Ubiquitin scan: SH3 p16'),
    ('uby_SH3', 'A', 'Ubiquitin scan: SH3_yeast p16'),
    ('ub_CUE', 'A', 'Ubiquitin scan: CUE p16'),
    ('uby_CUE', 'A', 'Ubiquitin scan: CUE_yeast p16'),
    ('ub_OTU', 'A', 'Ubiquitin scan: OTU p16'),
    ('uby_OTU', 'A', 'Ubiquitin scan: OTU_yeast p16'),
    ('ub_RPN13', 'A', 'Ubiquitin scan: RPN13 p16'),
    ('uby_RPN13', 'A', 'Ubiquitin scan: RPN13_yeast p16'),
    ('ub_UQcon', 'A', 'Ubiquitin scan: UQ_con p16'),
    ('uby_UQcon', 'A', 'Ubiquitin scan: UQ_con_yeast p16'),
]

def add_scores(output_dir, prediction_structure_scores_table, prediction_id_field, score_method_id):
    job_dict_path = os.path.join(os.path.join(output_dir, 'data'), 'job_dict.pickle')
    
    with open(job_dict_path, 'r') as f:
        job_dict = pickle.load(f)

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    db_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    DDGdb = db_api.DDG_db

    prediction_ids_and_structs_score_count = {}
    for row in DDGdb.execute_select("SELECT %s, ScoreType, StructureID FROM %s WHERE ScoreType IN ('WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex')" % (prediction_id_field, prediction_structure_scores_table)):
        prediction_id = long(row[prediction_id_field])
        score_type = row['ScoreType']
        structure_id = int(row['StructureID'])
        if (prediction_id, structure_id) not in prediction_ids_and_structs_score_count:
            prediction_ids_and_structs_score_count[(prediction_id, structure_id)] = 0
        prediction_ids_and_structs_score_count[(prediction_id, structure_id)] += 1
    structs_with_all_scores = set()
    for prediction_id, structure_id in prediction_ids_and_structs_score_count:
        if prediction_ids_and_structs_score_count[(prediction_id, structure_id)] == 6:
            structs_with_all_scores.add( (prediction_id, structure_id) )

    available_db3_files = {}
    available_db3_files_set = set()
    for task_name in job_dict:
        prediction_id = long(task_name.split('_')[0])
        round_num = int(task_name.split('_')[1])
        struct_type = task_name.split('_')[2]
        if struct_type == 'wt':
            continue
        wt_task_name = '%d_%d_%s' % (prediction_id, round_num, 'wt')
        mut_task_name = '%d_%d_%s' % (prediction_id, round_num, 'mut')
        wt_task_dir = os.path.join(output_dir, wt_task_name)
        mut_task_dir = os.path.join(output_dir, mut_task_name)
        wt_db3_file = os.path.join(wt_task_dir, 'output.db3.gz')
        mut_db3_file = os.path.join(mut_task_dir, 'output.db3.gz')
        if os.path.isfile(wt_db3_file) and os.path.isfile(mut_db3_file):
            available_db3_files_set.add( (prediction_id, round_num) )
            available_db3_files[(prediction_id, round_num)] = (wt_db3_file, mut_db3_file)
                
    db3_files_to_process = available_db3_files_set.difference(structs_with_all_scores)
    print 'Found %d output db3 files to add to score database' % len(db3_files_to_process)
    r = Reporter('parsing output db3 files and saving scores in database', entries='db3 files')
    r.set_total_count( len(db3_files_to_process) )
    scores_dict = {}
    for prediction_id, round_num in db3_files_to_process:
        wt_output_db3, mut_output_db3 = available_db3_files[(prediction_id, round_num)]

        tmp_dir = tempfile.mkdtemp(prefix='unzip_db3_')
        new_wt_output_db3_path = os.path.join(tmp_dir, os.path.basename(wt_output_db3))
        shutil.copy(wt_output_db3, new_wt_output_db3_path)
        wt_output_db3 = unzip_file(new_wt_output_db3_path)

        new_mut_output_db3_path = os.path.join(tmp_dir, os.path.basename(mut_output_db3))
        shutil.copy(mut_output_db3, new_mut_output_db3_path)
        mut_output_db3 = unzip_file(new_mut_output_db3_path)

        wtl_score = db_api.add_scores_from_db3_file(wt_output_db3, 1, round_num, db_api.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeLPartner', score_method_id = score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field))
        wtr_score = db_api.add_scores_from_db3_file(wt_output_db3, 2, round_num, db_api.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeRPartner', score_method_id = score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field))
        wtc_score = db_api.add_scores_from_db3_file(wt_output_db3, 3, round_num, db_api.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeComplex', score_method_id = score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field))
        ml_score = db_api.add_scores_from_db3_file(mut_output_db3, 1, round_num, db_api.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantLPartner', score_method_id = score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field))
        mr_score = db_api.add_scores_from_db3_file(mut_output_db3, 2, round_num, db_api.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantRPartner', score_method_id = score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field))
        mc_score = db_api.add_scores_from_db3_file(mut_output_db3, 3, round_num, db_api.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantComplex', score_method_id = score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field))
        scores_list = [wtl_score, wtr_score, wtc_score, ml_score, mr_score, mc_score]

        shutil.rmtree(tmp_dir)
        db_api.store_scores(None, prediction_id, scores_list, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
        r.increment_report()
    r.done()

def setup():
    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    db_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    DDGdb = db_api.DDG_db
    score_method_id = 7 # rescore with interface weights

    for ubiquitin_chain in ubiquitin_chains:
        results = DDGdb.execute_select('''
            SELECT PredictionSet, PDBFileID, COUNT( Prediction.ID ) AS NumRecords
            FROM `Prediction`
            INNER JOIN Experiment ON ExperimentID=Experiment.ID
            WHERE PredictionSet=%s
            GROUP BY PredictionSet
            ''', parameters = (ubiquitin_chain[2],))
        assert(len(results) == 1)
        assert(1425 <= results[0]['NumRecords'] <= 1445)

        cases = []
        results = DDGdb.execute_select('''
            SELECT Prediction.ID AS PredictionID, ProtocolID, Prediction.ExperimentID, PredictionSet, PDBFileID,
            Status, Errors
            FROM `Prediction`
            INNER JOIN Experiment ON ExperimentID=Experiment.ID
            WHERE PredictionSet=%s''', parameters = (ubiquitin_chain[2],))
        native = True
        organism = 'Human'
        if ubiquitin_chain[0].find('uby_') != 0:
            organism = 'Yeast'
            native = False
        case_name = ubiquitin_chain[0].replace('ub_', '').replace('uby_', '')
        if case_name == '1ubq':
            continue
        if organism == 'Human':
            continue
        for r in results:
            mutations = DDGdb.execute_select('SELECT * FROM ExperimentMutation WHERE ExperimentID=%s', parameters=(r['ExperimentID'],))
            assert(len(mutations) == 1)
            mutation = mutations[0]
            assert(r['Status'] == 'done' and r['Errors'] == None)
            prediction_id = r['PredictionID']
            archive = '/kortemmelab/shared/DDG/jobs/{0}.zip'.format(prediction_id)
            assert(os.path.exists(archive))
            cases.append(dict(
                case_name = case_name,
                organism = organism,
                native = native,
                prediction_id = r['PredictionID'],
                protocol_id = r['ProtocolID'],
                experiment_id = r['ExperimentID'],
                prediction_set = r['PredictionSet'],
                pdb_file_id = r['PDBFileID'],
                zip_filepath = archive,
                chain = mutation['Chain'],
                mutant_aa = mutation['MutantAA'],
                wildtype_aa = mutation['WildTypeAA'],
                residue_id = PDB.ResidueID2String(mutation['ResidueID']),
            ))

        already_unzipped_dirs = {}
        for d in os.listdir(tmpdir_location):
            if d.startswith('analyze_old_ubq_'):
                prediction_id = long(d.split('_')[3])
                already_unzipped_dirs[prediction_id] = os.path.join(tmpdir_location, d)
        
        print '%d total cases for %s' % ( len(cases), ubiquitin_chain[2] )
        r = Reporter('processing cases', entries='cases')
        r.set_total_count( len(cases) )

        for case in cases:
            prediction_id = long(case['prediction_id'])
            if prediction_id in already_unzipped_dirs:
                unzip_dir = already_unzipped_dirs[prediction_id]
            else:
                unzip_dir = tempfile.mkdtemp(prefix='analyze_old_ubq_%d_' % prediction_id, dir=tmpdir_location)
            ddg_output_path = os.path.join(unzip_dir, str(case['prediction_id']))
            if not prediction_id in already_unzipped_dirs:
                subprocess.check_output(['unzip', case['zip_filepath'], '-d', unzip_dir])
            db_api.add_rescore_cluster_run(
                ddg_output_path,
                case['chain'],
                score_method_id,
                prediction_id,
            )
            r.increment_report()
        r.done()

    db_api.create_cluster_run_rescore_dir( os.path.join(tmpdir_location, 'cluster_run') )
            
    # ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    # score_method_id = 7 # rescore with interface weigths
    # ppi_api.extract_data(prediction_set_name, root_directory = job_dir, score_method_id = score_method_id)

def add_scores_for_ubq_complex_runs()
    output_dir = '/dbscratch/kyleb/tmp/cluster_run/151112-kyleb_rescore_ddg_monomer'
    prediction_structure_scores_table = 'PredictionStructureScore'
    prediction_id_field = 'PredictionID'
    score_method_id = 7 # rescore with interface weights
    add_scores(output_dir, prediction_structure_scores_table, prediction_id_field, score_method_id)
    
if __name__ == '__main__':
    add_scores_for_ubq_complex_runs()
