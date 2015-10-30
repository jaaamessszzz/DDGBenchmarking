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
from klab.cluster_template.write_run_file import process as write_run_file
from subprocess import Popen, PIPE, check_call
import signal

job_output_directory = 'job_output'

class LineReader:
    def __init__(self,fname):
        if fname.endswith('.gz'):
            if not os.path.isfile(fname):
                raise IOError(fname)
            self.f = Popen(['gunzip', '-c', fname], stdout=PIPE, stderr=PIPE)
            self.zipped=True
        else:
            self.f = open(fname,'r')
            self.zipped=False
    def readlines(self):
        if self.zipped:
            for line in self.f.stdout:
                yield line
        else:
            for line in self.f.readlines():
                yield line
    def close(self):
        if self.zipped:
            if self.f.poll() == None:
                os.kill(self.f.pid,signal.SIGHUP)
        else:
            self.f.close()
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        self.close()
    def __iter__(self):
        return self.readlines()

def make_cst_file(cst_file_path, rosetta_outfile):
    if not os.path.isfile(rosetta_outfile):
        raise Exception("Couldn't find Rosetta outfile: " + rosetta_outfile)
    '''This does the work of convert_to_cst_file.sh'''

    constraints = []
    with LineReader(rosetta_outfile) as outf:
        for line in outf:
            if line.startswith("c-alpha"):
                line = line.split()
                constraints.append("AtomPair CA %s CA %s HARMONIC %s %s" % (line[5], line[7], line[9], line[12]))
    constraints = '\n'.join(constraints)
    # write constraints out
    with open(cst_file_path, 'w') as f:
        f.write(constraints)

if __name__ == '__main__':
    # Change this to match previous run
    prediction_set_id = 'ddg_monomer_16_002'

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path)

    prediction_ids = ppi_api.get_prediction_ids(prediction_set_id)

    existing_job = False
    end_job_name  = '%s_%s' % (getpass.getuser(), prediction_set_id)
    for d in os.listdir(job_output_directory):
        if os.path.isdir(os.path.join(job_output_directory, d)) and end_job_name in d:
            print 'Found existing job:', d
            job_name = d
            existing_job = True
    assert( existing_job )

    output_dir = os.path.join(job_output_directory, job_name )

    settings['scriptname'] = prediction_set_id + '_run-2'
    settings['tasks_per_process'] = 1
    settings['numjobs'] = len(prediction_ids)
    settings['mem_free'] = '4.1G'
    settings['output_dir'] = output_dir
    settings['job_dict_name'] = 'job_dict-2.pickle'
    settings['rosetta_args_list'] = [
        '-ignore_unrecognized_res',
        '-in:file:fullatom',
        '-fa_max_dis 9.0',
        '-ddg::dump_pdbs true',
        '-ddg::suppress_checkpointing true',
        '-ddg:weight_file soft_rep_design',
        '-ddg::iterations 50',
        '-ddg::local_opt_only false',
        '-ddg::min_cst true',
        '-ddg::mean false',
        '-ddg::min true',
        '-ddg::sc_min_only false',
        '-ddg::ramp_repulsive true',
    ]

    job_dict = {}

    prediction_ids = sorted( ppi_api.get_prediction_ids(prediction_set_id) )
    output_data_dir = os.path.join(settings['output_dir'], 'data')

    for prediction_id in prediction_ids:
        # Check if job already ran
        prev_prediction_id = prediction_id
        prev_prediction_id_dir = os.path.join(output_dir, str(prediction_id))
        prediction_id = '%d-ddg' % prediction_id
        prediction_id_dir = os.path.join(output_dir, prediction_id)
        if existing_job:
            rosetta_output_file = os.path.join( prediction_id_dir, 'rosetta.out.gz' )
            if os.path.isfile(rosetta_output_file):
                print 'Skipping', prediction_id
                settings['numjobs'] = settings['numjobs'] - 1
                continue
            if os.path.isdir(prediction_id_dir):
                print 'Job directory %s already exists, deleting' % prediction_id_dir
                shutil.rmtree(prediction_id_dir)
            # else:
            #     print 'Creating new job directory %s' % prediction_id_dir

        app_name = 'ddg_monomer'
        if 'appname' not in settings:
            settings['appname'] = app_name
        else:
            assert( settings['appname'] == app_name )

        job_details = ppi_api.get_job_details(prev_prediction_id)
        file_tuples = [] # List of names, contents
        for file_info in job_details['Files']['Input']:
            file_tuples.append( (file_info['Filename'], file_info['Content']) )
        substitution_parameters = json.loads(job_details['JSONParameters'])
        job_data_dir = os.path.join(output_data_dir, str(prev_prediction_id))

        files_dict = {} # Maps name to filepath position
        for file_name, file_contents in file_tuples:
            new_file_location = os.path.join(job_data_dir, file_name)
            with open(new_file_location, 'w') as f:
                f.write(file_contents)
            files_dict[file_name] = os.path.relpath(new_file_location, settings['output_dir'])
            
        cst_pdb = os.path.join(prev_prediction_id_dir, 'min_cst_0.5.%s' % substitution_parameters['%%input_pdb%%'])
        cst_pdb = cst_pdb[:-4] # strip .pdb
        cst_pdb += '_0001.pdb.gz'
        cst_file_path = os.path.join(job_data_dir, 'constraints.cst')
        make_cst_file(cst_file_path, os.path.join(prev_prediction_id_dir, 'rosetta.out.gz'))
        rel_cst_path = os.path.relpath(cst_file_path, settings['output_dir'])

        argdict = {
            '-in:file:s' : os.path.relpath(cst_pdb, output_dir),
            '-ddg::mut_file' : files_dict[substitution_parameters['%%pathtomutfile%%']],
            '-constraints::cst_file' : rel_cst_path,
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
