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
import time

def process_ddg_monomer_directory( job_dir ):
    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    date, prediction_set_name = os.path.basename(job_dir).split('-')
    # Strip user name off prediction_set_name
    if getpass.getuser() in prediction_set_name:
        prediction_set_name = prediction_set_name.strip( getpass.getuser() + '_' )
    
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    score_method_id = 8
    prediction_set_credit = 'Kyle Barlow'

    # Test get_top_x
    prediction_ids = ppi_api.get_prediction_ids(prediction_set_name)
    
    t1 = time.time()
    ppi_api.get_analysis_dataframe(prediction_set_name,
            prediction_set_series_name = 'Analysis1', prediction_set_description = 'Analysis1', prediction_set_credit = prediction_set_credit,
            use_existing_benchmark_data = True,
            include_derived_mutations = False,
            use_single_reported_value = False,
            take_lowest = 3,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            report_analysis = True,
            silent = False,
            root_directory = None,
            score_method_id = score_method_id,
            expectn = None,
            allow_failures = True,
            extract_data_for_case_if_missing = False,
            )

    # todo: store credit in dataframe or store/read from database
    ppi_api.analyze([prediction_set_name], score_method_id,
            analysis_set_ids = [],
            prediction_set_series_names = {}, prediction_set_descriptions = {}, prediction_set_credits = {}, prediction_set_colors = {}, prediction_set_alphas = {},
            use_existing_benchmark_data = True, recreate_graphs = False,
            include_derived_mutations = False,
            use_single_reported_value = False,
            expectn = 50,
            take_lowest = 3,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            output_directory = '/home/kyleb/tmp_analysis',
            generate_plots = True,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            )
    print('Time', time.time() - t1)

    #benchmark_run.calculate_metrics('ZEMu')


    sys.exit(0)
    
if __name__ == '__main__':
    job_dir = sys.argv[1]
    assert( os.path.isdir( job_dir ) )
    process_ddg_monomer_directory( job_dir )
