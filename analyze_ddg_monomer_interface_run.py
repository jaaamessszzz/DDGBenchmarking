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
    #### Change this
    prediction_set_name = ''

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']
    
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    score_method_id = 7 # rescore with interface weights score method (this is probably best...need to benchmark more)


    #### Uncomment the line matching the behavior you want:
    
    # To rescore the ddG output files on the fly:
    # root_directory should point to where the output files can be found (either the zipped file directories in /kortemmelab/shared or another location
    ppi_api.extract_data(prediction_set_name, root_directory = '/kortemmelab/shared/DDG/ppijobs', score_method_id = score_method_id)#, max_prediction_ids_to_process = 40)
    
    # To setup cluster run for rescoring:
    # ppi_api.extract_data(prediction_set_name, root_directory = '/kortemmelab/shared/DDG/ppijobs', score_method_id = score_method_id, setup_cluster_run = True)
    
    # To load into database rescoring cluster run result:
    # process_cluster_run(ppi_api, prediction_set_name, score_method_id, settings)

def process_cluster_run(ppi_api, prediction_set_name, score_method_id, settings):
    #### Set this to be a list of all output directories containing rescored PDBs
    #### There will probably be multiple output directories created by the setup cluster rescore function, to split up files and prevent rsync from failing
    #### Example: output_dirs = ['/dbscratch/kyleb/tmp/cluster_run/DiPUBSComplexes1-%d' % x for x in xrange(1, 12)]
    output_dirs = []

    # Leave these as is if your run is a normal one with results to be saved in PredictionPPIID (hack for one PUBS case)
    # This should be default argument for function below...but is not yet
    prediction_structure_scores_table = 'PredictionPPIStructureScore'
    prediction_id_field = 'PredictionPPIID'
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    ppi_api.add_scores_from_cluster_rescore(output_dirs, prediction_structure_scores_table, prediction_id_field, score_method_id)
    
if __name__ == '__main__':
    job_dir = sys.argv[1]
    assert( os.path.isdir( job_dir ) )
    process_ddg_monomer_directory( job_dir )
