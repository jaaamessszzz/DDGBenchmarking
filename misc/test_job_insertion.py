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


if __name__ == '__main__':
    ppi_api = get_ppi_interface(read_file('ddgdb.pw'),
                                rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease',
                                rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database')

    #pprint.pprint(ppi_api.get_score_method_details())

    #details = ppi_api.get_prediction_set_case_details('ZEMu run 1')
    #print(len(details['Data']))


    stability_api = get_protein_stability_interface(read_file('ddgdb.pw'))
    #pprint.pprint(stability_api.get_prediction_scores(55808))

    score_method_id = ppi_api.get_score_method_id('global', method_authors = 'kellogg')
    print(stability_api.get_top_x_ddg(55808, score_method_id, expectn = 49))

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
