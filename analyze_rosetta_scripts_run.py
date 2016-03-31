import os, sys
import shutil

import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
from klab.cluster_template.write_run_file import process as write_run_file
from ddglib.ppi_api import get_interface_with_config_file
from ddglib.rosetta_scripts_ppi_api import get_interface as get_interface_factory
import importlib

def main():
    assert( len(sys.argv) >= 3 )
    import_module = sys.argv[1]
    if import_module.endswith('.py'):
        import_module = import_module[:-3]
    if '/' in import_module:
        import_module = import_module.replace('/', '.')
    cfg = importlib.import_module(import_module, package=None)

    assert( os.path.isdir(sys.argv[2]) )
    root_directory = sys.argv[2]

    prediction_set_name = cfg.prediction_set_id

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )

    ppi_api.extract_data(prediction_set_name, root_directory = root_directory)

if __name__ == '__main__':
    main()
