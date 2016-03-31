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
import importlib

def process_ddg_monomer_directory():
    assert( len(sys.argv) >= 3 )
    import_module = sys.argv[1]
    if import_module.endswith('.py'):
        import_module = import_module[:-3]
    if '/' in import_module:
        import_module = import_module.replace('/', '.')
    cfg = importlib.import_module(import_module, package=None)
    processing_option = sys.argv[2]
    possible_options = ['on-the-fly', 'setup-rescore-run', 'load-rescore-run']

    if len(sys.argv) > 3 and os.path.isdir(sys.argv[3]):
        root_directory = sys.argv[3]
    else:
        root_directory = None

    prediction_set_name = cfg.prediction_set_id

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']
    
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    score_method_ids = cfg.score_method_ids

    for score_method_id in score_method_ids:
        if len(score_method_ids) > 1:
            print 'Processing score_method_id: %d\n' % score_method_id
        if processing_option == possible_options[0]:
            # To rescore the ddG output files on the fly:
            # root_directory should point to where the output files can be found (either the zipped file directories in /kortemmelab/shared or another location
            ppi_api.extract_data(prediction_set_name, root_directory = root_directory, score_method_id = score_method_id)#, max_prediction_ids_to_process = 40)

        elif processing_option == possible_options[1]:
            # To setup cluster run for rescoring:
            ppi_api.extract_data(prediction_set_name, root_directory = root_directory, score_method_id = score_method_id, setup_cluster_run = True)

        elif processing_option == possible_options[2]:
            # To load into database rescoring cluster run result:
            process_cluster_run(ppi_api, prediction_set_name, score_method_id, settings, output_dirs = [root_directory] )

        else:
            print 'ERROR: Argument 2 must be processing option. Choices:', possible_options

def process_cluster_run(ppi_api, prediction_set_name, score_method_id, settings, output_dirs = []):
    #### Set this to be a list of all output directories containing rescored PDBs
    #### There will probably be multiple output directories created by the setup cluster rescore function, to split up files and prevent rsync from failing
    #### Example: output_dirs = ['/dbscratch/kyleb/tmp/cluster_run/DiPUBSComplexes1-%d' % x for x in xrange(1, 12)]

    # Leave these as is if your run is a normal one with results to be saved in PredictionPPIID (hack for one PUBS case)
    # This should be default argument for function below...but is not yet
    prediction_structure_scores_table = 'PredictionPPIStructureScore'
    prediction_id_field = 'PredictionPPIID'
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    ppi_api.add_scores_from_cluster_rescore(output_dirs, prediction_structure_scores_table, prediction_id_field, score_method_id)
    
if __name__ == '__main__':
    process_ddg_monomer_directory()
