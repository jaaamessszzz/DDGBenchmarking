import os, sys
import shutil

import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
from klab.cluster_template.write_run_file import process as write_run_file
from ddglib.ppi_api import get_interface_with_config_file
from ddglib.ddg_monomer_ppi_api import get_interface as get_interface_factory

def process_ddg_monomer_directory( job_dir ):
    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    date, prediction_set_name = os.path.basename(job_dir).split('-')
    # Strip user name off prediction_set_name
    if getpass.getuser() in prediction_set_name:
        prediction_set_name = prediction_set_name.strip( getpass.getuser() + '_' )
    
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    score_method_id = 7 # rescore with interface weigths
    ppi_api.extract_data(prediction_set_name, root_directory = job_dir, score_method_id = score_method_id)
    
if __name__ == '__main__':
    job_dir = sys.argv[1]
    assert( os.path.isdir( job_dir ) )
    process_ddg_monomer_directory( job_dir )
