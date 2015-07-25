import os, sys

# Add parent directory to path
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from ddglib.ppi_api import get_interface_with_config_file
import tools.cluster_template.parse_settings as parse_settings
import time
import getpass
from tools.cluster_template.write_run_file import process as write_run_file

job_output_directory = 'job_output'

if __name__ == '__main__':
    # Change these for each run
    prediction_set_id = 'minimize_pack_bound_and_unbound_3cycles-1'
    script_file = 'minimize_pack_bound_and_unbound_3cycles.xml'
    
    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path)
    ppi_api.add_prediction_set(prediction_set_id, halted = True, priority = 7, allow_existing_prediction_set = True)

    # Populate the prediction set with jobs from a (tagged subset of a) user dataset
    ppi_api.add_prediction_run(prediction_set_id, 'AllBindingAffinity', tagged_subset = 'ZEMu', extra_rosetta_command_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res', show_full_errors = True)

    prediction_ids = ppi_api.get_prediction_ids(prediction_set_id)
    
    # for prediction_id in prediction_ids:
    #     details = ppi_api.get_job_details(prediction_id)
    #     # ppi_api.get_chains_for_mutatagenesis(details['PPMutagenesisID'], pdb_file_id, pdb_set_number, complex_id = None)
    #     print details.keys()
    #     print details['PDBMutations']
    #     for x in details['Files']['Input']:
    #         for key in x.keys():
    #             print '', key, x[key]
    #         print x['FileRole'], x['Filetype'], x['Filename']
    #     sys.exit(0)

    ppi_api.add_development_protocol_command_lines(
        prediction_set_id, prediction_set_id, 'rosetta_scripts',
        '-parser:protocol %s -in:file:s %%input_pdb%% -parser:script_vars chainstomove=%%chainstomove%% pathtoresfile=%%pathtoresfile%% -parser:view -inout:dbms:mode sqlite3 -inout:dbms:database_name rosetta_output.db3' % script_file,
        rosetta_script_file = 'interface/' + script_file,
    )

    job_name = '%s-%s_%s' % (time.strftime("%y%m%d"), getpass.getuser(), prediction_set_id)
    output_dir = os.path.join(job_output_directory, job_name )
    output_data_dir = os.path.join(output_dir, 'data')
    
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)

    settings['scriptname'] = prediction_set_id + '_run'
    settings['tasks_per_process'] = 1
    settings['numjobs'] = '%d' % len(prediction_ids)
    settings['mem_free'] = '1.2G'
    settings['output_dir'] = output_dir
    settings['db_id'] = prediction_set_id
    settings['run_from_database'] = True

    write_run_file(settings, database_run = True)

    print 'Job files written to directory:', os.path.abspath(output_dir)
    
# Unnecessary but here is how to change the values of batch_size, priority
# ppi_api.alter_prediction_set_batch_size(prediction_set_id, 40)
# ppi_api.alter_prediction_set_priority(prediction_set_id, 5)

# This should be called before kicking off jobs (or set halted = False above)
#ppi_api.start_prediction_set(prediction_set_id)

# compile the python submission script
