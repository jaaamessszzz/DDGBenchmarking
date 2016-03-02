#!/bin/python

import os, sys
import shutil
from klab.Reporter import Reporter
import zipfile
import multiprocessing

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

def zip_dir(output_path, ddg_dir, output_files):
    ddg_dir_path = os.path.join(output_path, ddg_dir)
    prediction_id, ddg_str = ddg_dir.split('-')
    prediction_id = long(prediction_id)
    all_files = find_all_files(ddg_dir_path)
    min_dir_path = os.path.join(output_path, '%d' % prediction_id)
    if os.path.isdir(min_dir_path):
        min_dir_files = find_all_files(min_dir_path)
    else:
        min_dir_files = []
    zipfile_path = os.path.join(ppijobs, '%d.zip' % prediction_id)
    if os.path.isfile(zipfile_path):
        print zipfile_path, 'already exists'
        return
    with zipfile.ZipFile(zipfile_path, 'w') as job_zip:
        for f in all_files:
            f_path = os.path.join(ddg_dir_path, f)
            job_zip.write(f_path, arcname=f)
        for output_file in output_files[prediction_id]:
            f_path = os.path.join(output_path, output_file)
            job_zip.write(f_path, arcname=os.path.join('output_logs', output_file) )
        for f in min_dir_files:
            f_path = os.path.join(ddg_dir_path, os.path.join(min_dir_path, f))
            job_zip.write(f_path, arcname=os.path.join('%d-premin' % prediction_id, f))

def parse_output_files(output_path):
    output_files = [f for f in os.listdir(output_path) if '.o' in f]
    output_files.extend( [f for f in os.listdir(output_path) if '.e' in f] )
    mapping_dict = {}
    for fpath in output_files:
        with open(os.path.join(output_path, fpath), 'r') as f:
            for line in f:
                if line.startswith('Making new job output directory'):
                    folder = line.strip().split('/')[-1].strip()
                    if folder.endswith('-ddg'):
                        folder = folder[:-4]
                    pred_id = long(folder)
                    if pred_id not in mapping_dict:
                        mapping_dict[pred_id] = []
                    mapping_dict[pred_id].append( fpath )
    return mapping_dict

def zip_directories(output_path, use_multiprocessing = True):
    ddg_dirs = [f for f in os.listdir(output_path) if f.endswith('-ddg')]
    output_files = parse_output_files(output_path)
    r = Reporter('zipping dirs', entries = 'dirs')
    r.set_total_count( len(ddg_dirs) )
    if use_multiprocessing:
        p = multiprocessing.Pool( min(4, multiprocessing.cpu_count()) )
    for ddg_dir in ddg_dirs:
        if use_multiprocessing:
            p.apply_async(zip_dir, (output_path, ddg_dir, output_files,), callback = r.increment_report_callback)
        else:
            r.increment_report_callback( zip_dir(output_path, ddg_dir, output_files) )
    if use_multiprocessing:
        p.close()
        p.join()
    r.done()

if __name__ == '__main__':
    assert( os.path.isdir(sys.argv[1]) )
    zip_directories( os.path.abspath(sys.argv[1]) )
