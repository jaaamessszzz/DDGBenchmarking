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

import tools.colortext as colortext
from tools.fs.fsio import read_file
from ddglib.ppi_api import get_interface as get_ppi_interface
from ddglib.monomer_api import get_interface as get_protein_stability_interface


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

    for prediction_record in ppi_api.get_queued_jobs(prediction_set_id, order_by = 'Cost', order_order_asc = False, include_files = True, truncate_content = 30):
        pprint.pprint(prediction_record['Cost'])

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
