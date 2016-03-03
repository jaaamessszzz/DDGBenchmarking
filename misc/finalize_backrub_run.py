import os, sys
import shutil

# Add parent directory to path
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
import shutil
import importlib
import imp
from klab.cluster_template.write_run_file import process as write_run_file
import cPickle as pickle

if __name__ == '__main__':
    assert( len(sys.argv) > 1 )
    cfg = imp.load_source('cfg', sys.argv[1])
    assert( os.path.isdir(sys.argv[2]) ) # Output directory

    prediction_set_id = cfg.prediction_set_id

    output_dir = sys.argv[2]
    data_dir_path = os.path.join(output_dir, 'data')
    job_dict_path = os.path.join(data_dir_path, 'blank_job_dict.pickle')
    settings_path = os.path.join(data_dir_path, 'settings.pickle')

    with open(job_dict_path, 'r') as f:
        job_dict = pickle.load(f)
    with open(settings_path, 'r') as f:
        settings = pickle.load(f)

    settings['output_dir'] = output_dir
    settings['tasks_per_process'] = 1
    settings['mem_free'] = '1.6G'
    settings['rosetta_args_list'].extend( [
        '-nstruct 50',
        '-ignore_unrecognized_res',
        '-ignore_zero_occupancy false',
        '-ex1', '-ex2',
        '-extrachi_cutoff 0',
        '-out:prefix bkrb_',
        '-mute core.io.pdb.file_data',
        '-backrub:ntrials 50000',
        '-mc_kt %.1f' % cfg.backrub_temp,
    ] )
    settings['appname'] = 'backrub'

    print job_dict
    print settings
    write_run_file(settings, database_run = False, job_dict = job_dict)
    print 'Job files written to directory:', os.path.abspath(output_dir)
