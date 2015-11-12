import os, sys
import shutil
import pprint

import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import tempfile
import json
import re
import datetime
import subprocess
from klab.cluster_template.write_run_file import process as write_run_file
from ddglib.rescore_db_api import get_interface as get_interface_factory
from ddglib.ppi_api import get_interface_with_config_file

from klab import colortext
from klab.bio.pdb import PDB

from ddglib import ddgdbapi, db_api

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

def main():
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

        print '%d cases for %s' % ( len(cases), ubiquitin_chain[2] )
        continue
        for case in cases:
            unzip_dir = tempfile.mkdtemp(prefix='analyze_old_ubq_')
            ddg_output_path = os.path.join(unzip_dir, str(case['prediction_id']))
            subprocess.check_output(['unzip', case['zip_filepath'], '-d', unzip_dir])
            print unzip_dir
            start = datetime.datetime.now()
            scores = db_api.parse_prediction_scores(
                case['prediction_id'],
                ddg_output_path=ddg_output_path,
                score_method_id=score_method_id,
                chains_to_move=case['chain']
            )
            print datetime.datetime.now() - start
            # print scores
            sys.exit(0)
            shutil.rmtree(unzip_dir)
            
    # ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    # score_method_id = 7 # rescore with interface weigths
    # ppi_api.extract_data(prediction_set_name, root_directory = job_dir, score_method_id = score_method_id)
    
if __name__ == '__main__':
    main()
