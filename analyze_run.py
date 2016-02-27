import matplotlib
matplotlib.use('Agg')
import os, sys
import shutil
import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
from klab.cluster_template.write_run_file import process as write_run_file
from klab.benchmarking.analysis.ddg_binding_affinity_analysis import DBBenchmarkRun as BindingAffinityBenchmarkRun
from ddglib.ppi_api import get_interface_with_config_file
from ddglib.ddg_monomer_ppi_api import get_interface as get_interface_factory
import datetime
import importlib
import tempfile

def process_ddg_monomer_directory():
    assert( len(sys.argv) >= 1 )

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = settings['local_rosetta_installation_path'] + '/database', get_interface_factory = get_interface_factory )

    all_benchmark_runs = []
    all_analysis_sets = set()
    for import_module in sys.argv[1:]:
        if import_module.endswith('.py'):
            import_module = import_module[:-3]
        if '/' in import_module:
            import_module = import_module.replace('/', '.')
        cfg = importlib.import_module(import_module, package=None)

        prediction_set_name = cfg.prediction_set_id
        prediction_set_credit = cfg.prediction_set_credit
        expectn = cfg.expectn
        take_lowests = cfg.take_lowests
        score_method_ids = cfg.score_method_ids

        benchmark_runs, analysis_sets = ppi_api.analyze(
            [prediction_set_name],
            score_method_ids,
            analysis_set_ids = ['ZEMu'],
            prediction_set_series_names = {}, prediction_set_descriptions = {}, prediction_set_colors = {}, prediction_set_alphas = {},
            use_existing_benchmark_data = cfg.use_existing_benchmark_data,
            recreate_graphs = False,
            include_derived_mutations = False,
            use_single_reported_value = False,
            prediction_set_credits = {prediction_set_name : prediction_set_credit},
            expectn = expectn,
            take_lowests = take_lowests,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            output_directory_root = None,
            generate_plots = True,
            generate_matplotlib_plots = True,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            debug = False,
            allow_failures = cfg.allow_missing_case_failures,
            call_analysis = False,
        )
        all_benchmark_runs.extend( benchmark_runs )
        all_analysis_sets.update( set(analysis_sets) )

    output_directory_root = '/tmp/%s/%s' % (getpass.getuser(), 'multiple_analysis')
    print 'Outputting to directory root:', output_directory_root, '\n'
    BindingAffinityBenchmarkRun.analyze_multiple(
        all_benchmark_runs,
        analysis_sets = list(all_analysis_sets),
        analysis_directory = output_directory_root,
        use_multiprocessing = True,
    )

if __name__ == '__main__':
    process_ddg_monomer_directory()
