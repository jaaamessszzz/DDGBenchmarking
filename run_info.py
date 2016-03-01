import os, sys
import shutil
import pprint
import pandas
import tempfile
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
from klab.Reporter import Reporter

job_output_directory = 'job_output'

def info_for_prediction_ids(prediction_ids):
    d = {
        'dataset_id' : [],
        'pdb' : [],
        'mutations' : [],
        'short_name' : [],
        ## 'long_name' : [],
    }
    r = Reporter('fetching job_details for prediction_ids', entries = 'prediction_ids')
    r.set_total_count( len(prediction_ids) )
    for prediction_id in prediction_ids:
        job_details = ppi_api.get_job_details(prediction_id)
        # rosetta_to_pdb_mapping = {long(key) : value for key, value in json.loads( [f['Content'] for f in job_details['Files']['Input'] if f['FileRole'] == 'Rosetta residue->PDB residue map'][0] ).iteritems()}
        d['dataset_id'].append( job_details['UserPPDataSetExperimentID'] )
        d['pdb'].append( job_details['Structure']['PDBFileID'] )
        d['short_name'].append( '%s / %s' % (job_details['Complex']['LShortName'], job_details['Complex']['RShortName']) )
        ## d['long_name'].append( '%s / %s' % (job_details['Complex']['LName'], job_details['Complex']['RName']) )
        mutations = ''
        for i, mutation in enumerate(job_details['PDBMutations']):
            residue = mutation['ResidueID'].strip()
            wt_aa = mutation['WildTypeAA']
            mut_aa = mutation['MutantAA']
            if i + 1 < len(job_details['PDBMutations']):
                mutations += ', %s%s%s' % (wt_aa, residue, mut_aa)
            else:
                mutations += '%s%s%s' % (wt_aa, residue, mut_aa)
        d['mutations'].append( mutations )
        r.increment_report()

    r.done()
    return pandas.DataFrame(d, index = prediction_ids)

if __name__ == '__main__':
    assert( len(sys.argv) > 1 )
    cfg = importlib.import_module(sys.argv[1], package=None)

    prediction_set_id = cfg.prediction_set_id
    protocol_name = cfg.protocol_name

    suppress_warnings = True

    # This uses the version of Rosetta from your cluster template settings file
    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database')

    assert( ppi_api.prediction_set_exists(prediction_set_id) )

    prediction_ids = sorted(ppi_api.get_prediction_ids(prediction_set_id))
    df = info_for_prediction_ids(prediction_ids)
    output_directory = tempfile.mkdtemp( prefix = '%s-%s-prediction_set_info_' % (time.strftime("%y%m%d"), getpass.getuser()) )
    output_csv = os.path.join(output_directory, '%s.csv' % prediction_set_id)
    print output_csv
    df.to_csv( output_csv )
