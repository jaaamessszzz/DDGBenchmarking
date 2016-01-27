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
    prediction_set_credit = cfg.prediction_set_credit
    expectn = cfg.expectn
    take_lowests = cfg.take_lowests
    score_method_ids = cfg.score_method_ids

    for score_method_id in score_method_ids:
        if len(score_method_ids) > 1:
            print 'Processing score method ID: %d' % score_method_id
        for take_lowest in take_lowests:
            if len(take_lowests) > 1:
                print 'Processing take_lowest (TopX): %d' % take_lowest
            output_dir_name = '%s-%s-%s_n-%d_topx-%d_score_method-%d_analysis' % (time.strftime("%y%m%d"), getpass.getuser(), prediction_set_name, expectn, take_lowest, score_method_id)

            output_directory = os.path.join('/tmp/%s/%s' % (getpass.getuser(), prediction_set_name), output_dir_name)

            if os.path.isdir( output_directory ):
                print 'Removing old output directory'
                shutil.rmtree( output_directory )

            print 'Outputting to directory:', output_directory, '\n'

            score_method_details = ppi_api.get_score_method_details(score_method_id = score_method_id)

            # todo: store credit in dataframe or store/read from database
            ppi_api.analyze([prediction_set_name], score_method_id,
                    analysis_set_ids = ['ZEMu'],
                    prediction_set_series_names = {}, prediction_set_descriptions = {}, prediction_set_credits = {prediction_set_name : '%s Score method id: %s (%d)' % (cfg.prediction_set_credit, score_method_details['MethodName'], score_method_id)}, prediction_set_colors = {}, prediction_set_alphas = {},
                    use_existing_benchmark_data = cfg.use_existing_benchmark_data,
                    recreate_graphs = False,
                    include_derived_mutations = False,
                    use_single_reported_value = False,
                    expectn = expectn,
                    take_lowest = take_lowest,
                    burial_cutoff = 0.25,
                    stability_classication_experimental_cutoff = 1.0,
                    stability_classication_predicted_cutoff = 1.0,
                    output_directory = output_directory,
                    generate_plots = True,
                    generate_matplotlib_plots = True,
                    report_analysis = True,
                    silent = False,
                    root_directory = None, # where to find the prediction data on disk
                    debug = False,
                    )
    
if __name__ == '__main__':
    process_ddg_monomer_directory()
