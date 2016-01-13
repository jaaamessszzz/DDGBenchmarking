#!/bin/python

import os, sys
import shutil
from klab.Reporter import Reporter
import zipfile
import multiprocessing

output_path = '/home/kyleb/data/ddg-benchmark-runs/150927-kyleb_ddg_monomer_16_002'
ppijobs = '/kortemmelab/shared/DDG/ppijobs'

def is_int(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def job_dir_from_output_file(output_file):
    job_dir_line_start = 'Job dir:'
    with open(output_file, 'r') as f:
        for line in f:
            if line.startswith(job_dir_line_start):
                return line[len(job_dir_line_start):].strip()
    raise Exception()

def move_output_files():
    output_files = []
    for f in os.listdir(os.getcwd()):
        if f.endswith('-ddg'):
            output_files.append(f)
        elif is_int(f):
            output_files.append(f)
    # output_files = [os.path.join(output_path, f) for f in os.listdir(output_path) if '.o' in f]
    r = Reporter('moving output files', entries='files')
    r.set_total_count(len(output_files))
    for output_file in output_files:
        shutil.move(output_file, os.path.join(output_path, job_dir_from_output_file(output_file)))
        r.increment_report()
    r.done()

def move_min_dirs():
    prediction_ids = [d for d in os.listdir(output_path) if is_int(d)]
    for prediction_id in prediction_ids:
        min_dir = os.path.join(output_path, prediction_id)
        ddg_dir = os.path.join(output_path, prediction_id + '-ddg')
        assert( os.path.isdir(min_dir) )
        assert( os.path.isdir(ddg_dir) )
        shutil.move(min_dir, ddg_dir)

def find_all_files(search_dir, prepend_dir=None):
    return_list = []
    for f in os.listdir(search_dir):
        f_path = os.path.join(search_dir, f)
        if os.path.isfile(f_path):
            if prepend_dir:
                return_list.append( os.path.join(prepend_dir, f) )
            else:
                return_list.append(f)
        elif os.path.isdir(f_path):
            return_list.extend( find_all_files(f_path, prepend_dir=os.path.basename(f_path)) )
    return return_list

def zip_dir(ddg_dir):
    ddg_dir_path = os.path.join(output_path, ddg_dir)
    prediction_id, ddg_str = ddg_dir.split('-')
    all_files = find_all_files(ddg_dir_path)
    with zipfile.ZipFile(os.path.join(ppijobs, '%s.zip' % prediction_id), 'w') as job_zip:
        for f in all_files:
            f_path = os.path.join(ddg_dir_path, f)
            job_zip.write(f_path, arcname=f)

def zip_directories():
    ddg_dirs = [f for f in os.listdir(output_path) if f.endswith('-ddg')]
    r = Reporter('zipping dirs', entries = 'dirs')
    r.set_total_count( len(ddg_dirs) )
    p = multiprocessing.Pool()
    for ddg_dir in ddg_dirs:
        p.apply_async(zip_dir, (ddg_dir,), callback = r.increment_report_callback)
    p.close()
    p.join()
    r.done()

if __name__ == '__main__':
    zip_directories()
