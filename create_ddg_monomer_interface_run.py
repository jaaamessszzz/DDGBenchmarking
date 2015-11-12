import os, sys
import shutil

# Add parent directory to path
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from ddglib.ppi_api import get_interface_with_config_file
import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
import shutil
from klab.cluster_template.write_run_file import process as write_run_file

job_output_directory = 'job_output'

def write_stripped_pdb(new_file_location, file_contents):
    if isinstance(file_contents, basestring):
        file_contents = file_contents.split('\n')
    with open(new_file_location, 'w') as f:
        for line in file_contents:
            if line.startswith('ATOM'):
                f.write(line + '\n')
    
if __name__ == '__main__':
    # Change these for each run
    prediction_set_id = 'DiPUBS: Complexes #1'
    protocol_name = 'ddg_monomer_16_002'

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database')
    # ppi_api.add_prediction_set(prediction_set_id, halted = True, priority = 7, allow_existing_prediction_set = True)

    # Populate the prediction set with jobs from a (tagged subset of a) user dataset
    # ppi_api.add_prediction_run(prediction_set_id, 'AllBindingAffinity', tagged_subset = 'ZEMu', extra_rosetta_command_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res', show_full_errors = True)

    prediction_ids = ppi_api.get_prediction_ids(prediction_set_id)

    # debugging lines
    # for prediction_id in prediction_ids:
    #     details = ppi_api.get_job_details(prediction_id)

    #     ppi_api.get_chains_for_mutatagenesis(details['PPMutagenesisID'], pdb_file_id, pdb_set_number, complex_id = None)
    #     print details.keys()
    #     print details['PDBMutations']
    #     for x in details['Files']['Input']:
    #         for key in x.keys():
    #             print '', key, x[key]
    #         print x['FileRole'], x['Filetype'], x['Filename']
    #     sys.exit(0)

    existing_job = False
    end_job_name  = '%s_%s' % (getpass.getuser(), prediction_set_id)
    for d in os.listdir(job_output_directory):
        if os.path.isdir(os.path.join(job_output_directory, d)) and end_job_name in d:
            print 'Found existing job:', d
            job_name = d
            existing_job = True
    if not existing_job:
        job_name = '%s-%s' % (time.strftime("%y%m%d"), end_job_name)
    
        ppi_api.add_development_protocol_command_lines(
            prediction_set_id, protocol_name, 'minimize_with_cst', ''
        )
        # 2x because bugs
        ppi_api.add_development_protocol_command_lines(
            prediction_set_id, protocol_name, 'minimize_with_cst', ''
        )

    output_dir = os.path.join(job_output_directory, job_name )

    settings['scriptname'] = 'diubq-pubs' + '_run'
    settings['tasks_per_process'] = 5
    settings['numjobs'] = len(prediction_ids)
    settings['mem_free'] = '3.0G'
    settings['output_dir'] = output_dir
    settings['rosetta_args_list'] = [
        '-in:file:fullatom',
        '-ignore_unrecognized_res',
        '-fa_max_dis 9.0',
        '-ddg::harmonic_ca_tether 0.5',
        '-ddg::constraint_weight 1.0',
        '-ddg::out_pdb_prefix min_cst_0.5',
        '-ddg::sc_min_only false',
    ]

    # Now get run settings from database and save to pickle file
    job_dict = {}
    output_data_dir = os.path.join(settings['output_dir'], 'data')

    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)

    prediction_ids = sorted( ppi_api.get_prediction_ids(prediction_set_id) )

    for prediction_id in prediction_ids:
        # Check if job already ran
        prediction_id_dir = os.path.join(output_dir, str(prediction_id))
        if existing_job:
            if os.path.isdir( prediction_id_dir ):
                pdb_output_files = [x for x in os.listdir( prediction_id_dir ) if '.pdb' in x]
            else:
                pdb_output_files = []
            if len(pdb_output_files) >= 1:
                # print 'Skipping', prediction_id
                settings['numjobs'] = settings['numjobs'] - 1
                continue
            if os.path.isdir(prediction_id_dir):
                print 'Job directory %s already exists, deleting' % prediction_id_dir
                shutil.rmtree(prediction_id_dir)
            # else:
            #     print 'Creating new job directory %s' % prediction_id_dir

        app_name = 'minimize_with_cst'
        if 'appname' not in settings:
            settings['appname'] = app_name
        else:
            assert( settings['appname'] == app_name )

        job_details = ppi_api.get_job_details(prediction_id)
        file_tuples = [] # List of names, contents
        for file_info in job_details['Files']['Input']:
            file_tuples.append( (file_info['Filename'], file_info['Content']) )
        substitution_parameters = json.loads(job_details['JSONParameters'])

        job_data_dir = os.path.join(output_data_dir, str(prediction_id))
        if os.path.isdir(job_data_dir):
            shutil.rmtree(job_data_dir)
        os.makedirs(job_data_dir)

        files_dict = {} # Maps name to filepath position
        for file_name, file_contents in file_tuples:
            new_file_location = os.path.join(job_data_dir, file_name)
            if '.pdb' in file_name:
                write_stripped_pdb(new_file_location, file_contents)
            else:
                with open(new_file_location, 'w') as f:
                    f.write(file_contents)
            files_dict[file_name] = os.path.relpath(new_file_location, settings['output_dir'])

        # Figure out input fi
            
        argdict = {
            'input_file_list' : [files_dict[substitution_parameters['%%input_pdb%%']]],
        }
        job_dict[prediction_id] = argdict

    write_run_file(settings, database_run = False, job_dict = job_dict)

    print 'Job files written to directory:', os.path.abspath(output_dir)

# Unnecessary but here is how to change the values of batch_size, priority
# ppi_api.alter_prediction_set_batch_size(prediction_set_id, 40)
# ppi_api.alter_prediction_set_priority(prediction_set_id, 5)

# This should be called before kicking off jobs (or set halted = False above)
#ppi_api.start_prediction_set(prediction_set_id)

# compile the python submission script
