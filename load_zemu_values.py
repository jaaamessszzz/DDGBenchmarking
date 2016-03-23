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

    q = "SELECT UserPPAnalysisSet.ID, UserPPAnalysisSet.Subset, UserPPAnalysisSet.RecordNumber, UserPPAnalysisSet.PPMutagenesisID, PPMutagenesisMutation.RecordKey, UserPPDataSetExperiment.PDBFileID FROM UserPPAnalysisSet RIGHT JOIN PPMutagenesisMutation ON PPMutagenesisMutation.PPMutagenesisID=UserPPAnalysisSet.PPMutagenesisID LEFT JOIN UserPPDataSetExperiment ON UserPPDataSetExperiment.PPMutagenesisID=PPMutagenesisMutation.PPMutagenesisID WHERE Subset='ZEMu' ORDER BY UserPPAnalysisSet.ID"
    results = ppi_api.DDG_db.execute_select(q)
    db_data = {}
    for i, row in enumerate(results):
        mut_id = row['ID']
        pdb_id = row['PDBFileID']
        if pdb_id not in db_data:
            db_data[pdb_id] = {}
        if mut_id not in db_data[pdb_id]:
            db_data[pdb_id][mut_id] = set()
        record_key = row['RecordKey']
        wt_aa = record_key.split()[0]
        assert( len(wt_aa) == 1 )
        chain = record_key.split()[1][0]
        mut_aa = record_key.split()[1][-1]
        resnum = record_key.split()[1][1:-1]
        data_key = '%s%s%s%s' % (chain, wt_aa, resnum, mut_aa)
        db_data[pdb_id][mut_id].add( data_key )

    db_frozen_data = {}
    for pdb_id in db_data:
        db_frozen_data[pdb_id] = {}
        for k, v in db_data[pdb_id].iteritems():
            db_frozen_data[pdb_id][ frozenset(v) ] = k

    zemu_data_file = 'misc/annotated_main_zemu_dataset.csv'
    with open(zemu_data_file, 'r') as f:
        for line in f.readlines()[1:]:
            data = line.strip().split('\t')
            zemu_id = long(data[0])
            pdb_id = data[1]
            mutations = frozenset([x.strip() for x in data[2].split(',')])
            if pdb_id not in db_frozen_data:
                print 'missing:', pdb_id
            elif mutations not in db_frozen_data[pdb_id]:
                print 'missing:', mutations
            else:
                db_id = db_frozen_data[pdb_id][mutations]
                print 'match:', db_id
    sys.exit(0)

    prediction_ids = sorted(ppi_api.get_prediction_ids(prediction_set_id))
    for prediction_id in prediction_ids:
        pass
