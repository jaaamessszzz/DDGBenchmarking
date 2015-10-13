import sys
import traceback
import datetime
import string
import os
import pickle
import numpy
import json
import pprint

if __name__ == "__main__":
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, "../updatedb")
    sys.path.insert(0, '/home/oconchus/dev/')

import tools.colortext as colortext
from tools.fs.fsio import read_file, write_file
from ddglib.ppi_api import get_interface as get_ppi_interface
from ddglib.ddg_monomer_ppi_api import get_interface as get_kyles_ppi_interface
from ddglib.monomer_api import get_interface as get_protein_stability_interface

def export_datasets():
    ppi_api = get_ppi_interface(read_file('ddgdb.pw'),
                                rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease',
                                rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database')
    dataset_set_ids = ['ZEMu_10.1002/prot.24634', 'SKEMPI_2012/10/03', 'BeAtMuSiC_10.1093/nar/gkt450']
    for dataset_set_id in dataset_set_ids:
        dataset_set_short_name = dataset_set_id.split('_')[0]
        colortext.warning(dataset_set_short_name)
        dataset_set_json = ppi_api.export_dataset_to_json(dataset_set_id)
        dataset_set_tsv = ppi_api.export_dataset_to_csv(dataset_set_id)
        write_file('%s.export.tsv' % dataset_set_short_name, dataset_set_tsv)
        write_file('%s.export.json' % dataset_set_short_name, dataset_set_json)


def add_dummy_data(ppi_api):
    '''A function to add random data to the database to test the analysis API.'''

    import random
    score_method_id = ppi_api.get_score_method_id('interface', method_authors = 'kyle')
    prediction_records = ppi_api.DDG_db.execute_select('SELECT * FROM PredictionPPI WHERE PredictionSet="ZEMu run 1"')
    for prediction_record in prediction_records:
        UserPPDataSetExperimentID = prediction_record['UserPPDataSetExperimentID']
        analysis_records = ppi_api.get_experimental_ddgs_by_analysis_set(user_dataset_experiment_id = UserPPDataSetExperimentID)

        experimental_value = analysis_records[UserPPDataSetExperimentID]['ZEMu']['MeanDDG']

        # Create 50 random predicted values
        jitter = random.uniform(-1, 1)
        adj_experimental_value = experimental_value + jitter
        for i in range(1, 51):
            predicted_mutant_value = adj_experimental_value + random.uniform(-0.3, 0.3)

            ppi_api.DDG_db.insertDictIfNew('PredictionPPIStructureScore', dict(
                PredictionPPIID = prediction_record['ID'],
                ScoreMethodID = score_method_id,
                ScoreType = 'DDG',
                StructureID = i,
                total = predicted_mutant_value
            ), ['PredictionPPIID', 'ScoreMethodID', 'ScoreType', 'StructureID'])
        sys.stdout.write('.'); sys.stdout.flush()


if __name__ == '__main__':
    ppi_api = get_ppi_interface(read_file('ddgdb.pw'),
                                rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease',
                                rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database')

    kyles_ppi_api = get_kyles_ppi_interface(read_file('ddgdb.pw'),
                                rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease',
                                rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database')


    #pprint.pprint(ppi_api.get_score_method_details())
    #details = ppi_api.get_prediction_set_case_details('ZEMu run 1')
    #print(len(details['Data']))

    prediction_set_id = 'ZEMu run 1'
    score_method_id = ppi_api.get_score_method_id('interface', method_authors = 'kyle', method_type = 'global')
    prediction_set_credit = "Shane O'Connor"

    #prediction_set_id = 'ddg_monomer_16_002'
    #score_method_id = kyles_ppi_api.get_score_method_id('Rescore-Interface', method_authors = 'kyle')
    #prediction_set_credit = 'Kyle Barlow'

    import time
    t1 = time.time()
    ppi_api.get_analysis_dataframe(prediction_set_id,
            prediction_set_series_name = 'My test', prediction_set_description = 'My test description', prediction_set_credit = prediction_set_credit,
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
            allow_failures = False,
            )

    # todo: store credit in dataframe or store/read from database
    ppi_api.analyze([prediction_set_id], score_method_id,
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
            output_directory = '/tmp/analysis',
            generate_plots = True,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            )
    print('Time', time.time() - t1)

    #benchmark_run.calculate_metrics('ZEMu')


    sys.exit(0)
    s1 = ppi_api.get_score_dict(prediction_id = 1265, score_method_id = '4', score_type = 'WildTypeLPartner', structure_id = '23')
    s2 = ppi_api.get_score_dict(prediction_id = 1265, score_method_id = '4', score_type = 'WildTypeRPartner', structure_id = '23')
    s3 = ppi_api.get_score_dict(prediction_id = 1265, score_method_id = '4', score_type = 'WildTypeComplex', structure_id = '24')
    print('r')
    ppi_api.store_scores('ZEMu run 1', 1265, [s1, s2, s3])
    a='''
    {1L: {1L: {'Mutant': {'DDG': None,
                      'dslf_ca_dih': None,
                      'dslf_cs_ang': None,
                      'dslf_fa13': 0.0,
                      ...
     4L: {'None': {'DDG': {'DDG': 0.887,
                       'dslf_ca_dih': None,
                       'dslf_cs_ang': None,
                       'dslf_fa13': None,

    '''


    ppi_api.get_analysis_dataframe('ZEMu run 1',
            prediction_set_series_name = 'My test', prediction_set_description = 'My test description', prediction_set_credit = 'Shane',
            use_existing_benchmark_data = False, # todo
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
            allow_failures = False,
            )

    print('s')

    #print(ppi_api.DDG_db.FieldNames.__dict__('PredictionPPIStructureScore'))#ppi_api.])


    #ScoreMethodID - 1,2,3,4
    #ScoreType - 'DDG', 'WildType', 'Mutant', 'WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex' VARCHAR(32)
    #StructureID -  nstruct 1-50
    #DDG - score term


sys.exit(0)

if __name__ == '__main__':
    ppi_api = get_ppi_interface(read_file('ddgdb.pw'),
                                rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease',
                                rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database')
    prediction_set_id = 'ZEMu run multiprocessing test'
    ppi_api.add_prediction_set(prediction_set_id, halted = True, priority = 7, batch_size = 41, allow_existing_prediction_set = True)
    ppi_api.alter_prediction_set_batch_size(prediction_set_id, 40)
    ppi_api.alter_prediction_set_priority(prediction_set_id, 5)
    ppi_api.start_prediction_set(prediction_set_id)
    #ppi_api.add_prediction_run(prediction_set_id, 'AllBindingAffinity', tagged_subset = 'ZEMu', extra_rosetta_command_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res', test_run_first = False, test_only = False, show_full_errors = True)
    ppi_api.add_prediction_run_mp(prediction_set_id, 'AllBindingAffinity', tagged_subset = 'ZEMu', extra_rosetta_command_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res', show_full_errors = True)
    sys.exit(0)

if __name__ == '__main__':

    # Set up the stability API
    stability_api = get_protein_stability_interface(read_file('ddgdb.pw'))

    # Set up the PPI API and get direct access to the database interface objects
    ppi_api = get_ppi_interface(read_file('ddgdb.pw'),
                                rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease',
                                rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database')
    ddg_db = ppi_api.DDG_db
    ddg_db_utf = ppi_api.DDG_db_utf

    # Informational_job tests
    ppi_api.get_prediction_set_details('RosCon2013_P16_score12prime')
    stability_api.get_prediction_ids('RosCon2013_P16_score12prime')

    # Print API help
    #ppi_api.help()

    # Create the prediction set
    prediction_set_id = 'ZEMu run 1'
    ppi_api.add_prediction_set(prediction_set_id, halted = True, priority = 7, batch_size = 41, allow_existing_prediction_set = True)
    ppi_api.alter_prediction_set_batch_size(prediction_set_id, 40)
    ppi_api.alter_prediction_set_priority(prediction_set_id, 5)
    ppi_api.start_prediction_set(prediction_set_id)

    # Set up the job specific details
    # This needs to be done since there are no associated protocols yet; Kyle is in charge of this part.
    #for prediction_id in ppi_api.get_prediction_ids(prediction_set_id):

    x = 0
    printers = [colortext.pgreen, colortext.pyellow]
    for job_details in ppi_api.get_queued_jobs(prediction_set_id, order_by = 'Cost', order_order_asc = False, include_files = True, truncate_content = 30):
        if x >= 3:
            break
        printers[x % 2](pprint.pformat(job_details))
        x += 1

    sys.exit(0)
    while True:
        job_details = ppi_api.get_job(prediction_set_id, order_by = 'Cost', order_order_asc = False, include_files = True)
        colortext.pyellow(pprint.pformat(job_details))
        break

    #    print(prediction_id)
    #    #... get job details, create command lines
    #    #ppi_api.add_job_command_lines(prediction_id, cmd_lines)

    sys.exit(0)

    # Print the available user datasets
    #pprint.pprint(ppi_api.get_defined_user_datasets())

    # Populate the prediction set with jobs from a (tagged subset of a) user dataset
    ppi_api.add_prediction_run(prediction_set_id, 'AllBindingAffinity', tagged_subset = 'ZEMu', extra_rosetta_command_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res', test_run_first = True, show_full_errors = True)

    # Set up the job specific details
    # This needs to be done since there are no associated protocols yet; Kyle is in charge of this part.
    for prediction_id in ppi_api.get_prediction_ids(prediction_set_id):
        print(prediction_id)
        #... get job details, create command lines
        #ppi_api.add_job_command_lines(prediction_id, cmd_lines)


    sys.exit(0)
    a='''
    ppi_api.start_prediction_set(prediction_set_id)

    # Run the jobs
    while True:
        j = ppi_api.get_job(prediction_set)
        if j:
            ...create the jobs files on disk, group into an array job
        else:
            break
    #submit the job to the cluster, calling start_job(prediction_id) when the jobs have been submitted to SGE

    Inserting a job

    '''
