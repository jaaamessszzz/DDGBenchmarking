import os, sys
import shutil

if __name__ == '__main__': # Hack for Shane: delete at will!
    sys.path.insert(0, '../klab')

import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
from klab.cluster_template.write_run_file import process as write_run_file
from ddglib.ppi_api import get_interface_with_config_file
from ddglib.ddg_monomer_ppi_api import get_interface as get_interface_factory
import datetime
import importlib
import tempfile


def process_ddg_monomer_directory():
    assert( len(sys.argv) >= 1 )
    import_module = sys.argv[1]
    if import_module.endswith('.py'):
        import_module = import_module[:-3]
    if '/' in import_module:
        import_module = import_module.replace('/', '.')
    cfg = importlib.import_module(import_module, package=None)

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    prediction_set_name = cfg.prediction_set_id
    
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = settings['local_rosetta_installation_path'] + '/database', get_interface_factory = get_interface_factory )
    score_method_id = cfg.score_method_id
    prediction_set_credit = cfg.prediction_set_credit

    # Test get_top_x
    prediction_ids = ppi_api.get_prediction_ids(prediction_set_name)
    
    t1 = datetime.datetime.now()
    # ppi_api.get_analysis_dataframe(prediction_set_name,
    #         prediction_set_series_name = 'DefaultAnalysis', prediction_set_description = 'DefaultAnalysis', prediction_set_credit = prediction_set_credit,
    #         use_existing_benchmark_data = True,
    #         include_derived_mutations = False,
    #         use_single_reported_value = False,
    #         take_lowest = 3,
    #         burial_cutoff = 0.25,
    #         stability_classication_experimental_cutoff = 1.0,
    #         stability_classication_predicted_cutoff = 1.0,
    #         report_analysis = True,
    #         silent = False,
    #         root_directory = None,
    #         score_method_id = score_method_id,
    #         expectn = None,
    #         allow_failures = True,
    #         extract_data_for_case_if_missing = False,
    #         )

    # todo: store credit in dataframe or store/read from database
    ppi_api.analyze([prediction_set_name], score_method_id,
            analysis_set_ids = ['ZEMu'],
            prediction_set_series_names = {}, prediction_set_descriptions = {}, prediction_set_credits = {prediction_set_name : 'kyleb'}, prediction_set_colors = {}, prediction_set_alphas = {},
            use_existing_benchmark_data = False, recreate_graphs = False,
            include_derived_mutations = False,
            use_single_reported_value = False,
            expectn = 45,
            take_lowest = 3,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            output_directory = tempfile.mkdtemp(prefix='%s-%s-%s-analysis_' % (time.strftime("%y%m%d"), getpass.getuser(), prediction_set_name) ),
            generate_plots = True,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            debug = False,
            )
    print('Time', datetime.datetime.now() - t1)

    #benchmark_run.calculate_metrics('ZEMu')


    sys.exit(0)
    
if __name__ == '__main__':
    process_ddg_monomer_directory()
