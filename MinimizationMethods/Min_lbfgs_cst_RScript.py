#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-1240
#$ -l arch=linux-x64
#$ -l mem_free=1.1G
#$ -l netapp=1G,scratch=1G

import socket
import sys

print "Python version:", sys.version
print "Hostname:", socket.gethostname()

from datetime import *
import os
import subprocess
import time
import sys
import shutil
import inspect
import gzip
import tempfile
import re
import json

def roundTime(dt=None, roundTo=1):
    """
    Round a datetime object to any time period (in seconds)
    dt : datetime.datetime object, default now.
    roundTo : Closest number of seconds to round to, default 1 second.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
    http://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python/10854034#10854034
    """
    if dt == None : dt = datetime.now()
    seconds = total_seconds(dt - dt.min)
    # // is a floor division, not a comment on following line:
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + timedelta(0,rounding-seconds,-dt.microsecond)

def total_seconds(td):
    '''
    Included in python 2.7 but here for backwards-compatibility for old Python versions (like on QB3 cluster)
    '''
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6

time_start = roundTime()

sge_task_id=0
if os.environ.has_key("SGE_TASK_ID"):
    sge_task_id = long(os.environ["SGE_TASK_ID"])

job_id=0
if os.environ.has_key("JOB_ID"):
    job_id=long(os.environ["JOB_ID"])

print 'Starting time:', time_start
print 'Job id:', job_id
print 'Task id:', sge_task_id
print os.getcwd()
# Run task (example code include below)

#Parses dataset .json file and outputs chain to move and input PDB file directory
def json_parser(datadir):
    jsonload = open("%sdata/blank_job_dict.json" %datadir)
    jsonfile = json.load(jsonload)
    key = sorted(jsonfile.keys())[sge_task_id-1]
    inputdir = jsonfile[key]['input_file_list'][0]
    return inputdir

def minimize(inputdir, outputdir):
    #Removes PDB file from path, saves in variable filenum_dir
    inputdir_parse = re.sub("/",' ', str(inputdir))
    data, filenum, pdbtemp = inputdir_parse.split()
    filenum_dir = os.path.join (data, filenum)
    PDBID = pdbtemp[:-4]
    predIDoutdir = os.path.join (outputdir, filenum)
    
    #Makes a folder for data dumping
    print 'Making directory %s%s...' %(outputdir, filenum)
    try:
        os.makedirs(predIDoutdir)
    except:
        '%s already exists!' % predIDoutdir
        
    #Assigns function output to variables for bash input (resfile_relpath)
    pdb_relpath = os.path.relpath('/netapp/home/james.lucas/160322-james-backrub-rscript-full/%s' %inputdir, predIDoutdir)
    xml_relpath = os.path.relpath('/netapp/home/james.lucas/DDGBenchmarks_Test/MinimizationMethods/', os.getcwd())
    print os.getcwd()
    
    arg = ['/netapp/home/james.lucas/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease',
           '-s',
           pdb_relpath,
           '-parser:protocol',
           '/netapp/home/james.lucas/DDGBenchmarks_Test/MinimizationMethods/Minimization_lbfgs_armijo_nonmonotone_cst.xml',
           '-ignore_unrecognized_res',
           '-nstruct 10'
          ]
    
    print 'Working on: %s %s' %(filenum, PDBID)
    print ' '.join(arg)
    
    outfile_path = os.path.join(predIDoutdir, 'rosetta.out')
    rosetta_outfile = open(outfile_path, 'w')
    
    print 'Running RosettaScript...'
    try:
        rosetta_process = subprocess.Popen(arg, stdout=rosetta_outfile, cwd=predIDoutdir)
    except OSError:
        print arg, rosetta_outfile, predIDoutdir
        raise
    return_code = rosetta_process.wait()
    
    print 'Task return code:', return_code, '\n'
    rosetta_outfile.close()    
    
    return filenum

datadir = '/netapp/home/james.lucas/160322-james-backrub-rscript-full/'
inputdir = json_parser(datadir)
outputdir = '/netapp/home/james.lucas/MinimizationMethods/output-Min_lbfgs_cst/'
filenum = minimize(inputdir, outputdir)
#End Pasted Stuff

time_end = roundTime()
print 'Ending time:', time_end
print "Elapsed time:", time_end-time_start

# Calculate RAM usage                                                             
qstat_p = subprocess.Popen(['/usr/local/sge/bin/linux-x64/qstat', '-j', '%d' % job_id],
                           stdout=subprocess.PIPE)
out, err = qstat_p.communicate()

for line in out.split(os.linesep):
    m = re.match('(?:usage\s+%d[:]\s+.*?)(?:maxvmem[=])(\d+[.]\d+)([a-zA-Z]+)(?:.*?)' % sge_task_id, line)
    if m:
        ram_usage = float(m.group(1))
        ram_usage_type = m.group(2)
        print 'Max virtual memory usage: %.1f%s' % (ram_usage, ram_usage_type)

#moves output and error files to output directory
#error_out = 'Min_lbfgs_cst_RScript.py.e' + str(job_id) + '.' + str(sge_task_id)
#output_out = 'Min_lbfgs_cst_RScript.py.o' + str(job_id) + '.' + str(sge_task_id)

#shutil.move(error_out , os.path.join(outputdir, filenum))
#shutil.move(output_out , os.path.join(outputdir, filenum))