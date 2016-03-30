import os, sys
import shutil
import pprint
import json

# Add parent directory to path
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from ddglib.ppi_api import get_interface_with_config_file
import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
import shutil
import importlib
from klab.cluster_template.write_run_file import process as write_run_file
from klab import colortext
from klab.fs.fsio import write_file
import cPickle as pickle

job_output_directory = 'job_output'


if __name__ == '__main__':
    assert( len(sys.argv) == 1 )
    cfg = importlib.import_module('run_config.zemu-values', package=None)

    prediction_set_id = cfg.prediction_set_id

    suppress_warnings = True

    # This uses the version of Rosetta from your cluster template settings file
    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database')

    zemu_data = {}
    zemu_data_file = 'misc/annotated_main_zemu_dataset.csv'
    with open(zemu_data_file, 'r') as f:
        for line in f.readlines()[1:]:
            data = line.strip().split('\t')
            zemu_id = long(data[0])
            pdb_id = data[1]
            zemu_data[zemu_id] = float(data[5])

    q = "SELECT PredictionPPI.ID, PredictionPPI.UserPPDataSetExperimentID, RecordNumber, Status, UserPPAnalysisSet.Subset FROM ddG.PredictionPPI INNER JOIN UserPPAnalysisSet ON UserPPAnalysisSet.UserPPDataSetExperimentID=PredictionPPI.UserPPDataSetExperimentID WHERE PredictionSet='zemu-values' AND UserPPAnalysisSet.Subset='ZEMu'"
    results = ppi_api.DDG_db.execute_select(q)
    db_data = {}
    for i, row in enumerate(results):
        ppi_id = row['ID']
        record_num = row['RecordNumber']
        score = ppi_api.get_score_dict(prediction_id = ppi_id, score_method_id = cfg.score_method_ids[0], score_type = 'DDG', structure_id = 1)
        score['total'] = zemu_data[record_num]
        ppi_api.store_scores( cfg.prediction_set_id, ppi_id, [score] )
