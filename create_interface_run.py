import os, sys

# Add parent directory to path
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from ddglib.ppi_api import get_interface_with_config_file
import tools.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
from tools.cluster_template.write_run_file import process as write_run_file

job_output_directory = 'job_output'
run_from_database = False # Controls if each cluster node attempts to get its info directly from the DB

if __name__ == '__main__':
    # Change these for each run
    prediction_set_id = 'pack_bound_and_unbound_3cycles-4'
    script_file = 'pack_bound_and_unbound.xml'

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
        '-parser:protocol ' + str(script_file) + ' -in:file:s %%input_pdb%% -parser:script_vars chainstomove=%%chainstomove%% pathtoresfile=%%pathtoresfile%% -parser:view -inout:dbms:mode sqlite3 -inout:dbms:database_name rosetta_output.db3',
        rosetta_script_file = 'interface/' + script_file,
    )
    # 2x because bugs
    ppi_api.add_development_protocol_command_lines(
        prediction_set_id, prediction_set_id, 'rosetta_scripts',
        '-parser:protocol ' + str(script_file) + ' -in:file:s %%input_pdb%% -parser:script_vars chainstomove=%%chainstomove%% pathtoresfile=%%pathtoresfile%% -parser:view -inout:dbms:mode sqlite3 -inout:dbms:database_name rosetta_output.db3',
        rosetta_script_file = 'interface/' + script_file,
    )

    job_name = '%s-%s_%s' % (time.strftime("%y%m%d"), getpass.getuser(), prediction_set_id)
    output_dir = os.path.join(job_output_directory, job_name )

    settings['scriptname'] = prediction_set_id + '_run'
    settings['tasks_per_process'] = 1
    settings['numjobs'] = '%d' % len(prediction_ids)
    settings['mem_free'] = '1.2G'
    settings['output_dir'] = output_dir
    settings['db_id'] = prediction_set_id

    if not run_from_database:
        # Now get run settings from database and save to pickle file
        job_dict = {}
        output_data_dir = os.path.join(settings['output_dir'], 'data')

        if not os.path.isdir(output_data_dir):
            os.makedirs(output_data_dir)

        prediction_ids = sorted( ppi_api.get_prediction_ids(prediction_set_id) )

        for task_id in xrange(0, int(settings['numjobs'])):
            prediction_id = prediction_ids[task_id]
            job_details = ppi_api.get_job_details(prediction_id)
            if not job_details['DevelopmentProtocolID']:
                raise Exception("Missing DevelopmentProtocolID")
            development_protocol = ppi_api.get_development_protocol(job_details['DevelopmentProtocolID'])
            app_name = development_protocol['Application']
            if 'appname' not in settings:
                settings['appname'] = app_name
            else:
                assert( settings['appname'] == app_name )
            flags_list = development_protocol['TemplateCommandLine'].strip().split()
            file_tuples = [] # List of names, contents
            for file_info in job_details['Files']['Input']:
                file_tuples.append( (file_info['Filename'], file_info['Content']) )
            substitution_parameters = json.loads(job_details['JSONParameters'])
            extra_parameters = job_details['ExtraParameters']
            job_data_dir = os.path.join(output_data_dir, str(prediction_id))
            if not os.path.isdir(job_data_dir):
                os.makedirs(job_data_dir)

            # Add extra parameters to flags_list
            for extra_parameter in extra_parameters.strip().split():
                flags_list.append(extra_parameter)

            files_dict = {} # Maps name to filepath position
            for file_name, file_contents in file_tuples:
                new_file_location = os.path.join(job_data_dir, file_name)
                with open(new_file_location, 'w') as f:
                    f.write(file_contents)
                files_dict[file_name] = os.path.relpath(new_file_location, settings['output_dir'])

            arglist = []
            argdict = {}
            parsing_scriptvars = False
            scriptvar_list = []
            for flag in flags_list:
                matches = re.findall('%%.+%%', flag)
                for match_str in matches:
                    if match_str in substitution_parameters:
                        if substitution_parameters[match_str] in files_dict:
                            flag = flag.replace(match_str, files_dict[substitution_parameters[match_str]])
                        else:
                            flag = flag.replace(match_str, substitution_parameters[match_str])

                if parsing_scriptvars:
                    if '=' in flag:
                        scriptvar_list.append(flag)
                    else:
                        parsing_scriptvars = False

                if flag == '-parser:script_vars':
                    parsing_scriptvars = True

                if not parsing_scriptvars:
                    # Check if argument is a file
                    if flag in files_dict:
                        last_arg = arglist.pop()
                        file_flag = files_dict[flag]
                        argdict[last_arg] = files_dict[flag]
                    else:
                        arglist.append(flag)

            if len(scriptvar_list) > 0:
                argdict['-parser:script_vars'] = scriptvar_list

            if len(arglist) > 0:
                argdict['FLAGLIST'] = arglist

            job_dict[prediction_id] = argdict
    else:
         job_dict = None

    write_run_file(settings, database_run = run_from_database, job_dict = job_dict)

    print 'Job files written to directory:', os.path.abspath(output_dir)

# Unnecessary but here is how to change the values of batch_size, priority
# ppi_api.alter_prediction_set_batch_size(prediction_set_id, 40)
# ppi_api.alter_prediction_set_priority(prediction_set_id, 5)

# This should be called before kicking off jobs (or set halted = False above)
#ppi_api.start_prediction_set(prediction_set_id)

# compile the python submission script
