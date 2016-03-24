#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-850
#$ -l arch=linux-x64
#$ -l mem_free=8G
#$ -l netapp=2G,scratch=1G

# Make sure you set task number above to be correct!!!

import socket
import sys

print "Python version:", sys.version
print "Hostname:", socket.gethostname()

from Bio.PDB import *
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
    
import json
import re

def read_mutations_resfile(filenum_dir):
    resfile = os.path.join(filenum_dir, 'mutations.resfile')
    mutations = []
    with open(resfile, 'r') as f:
        post_start = False
        for line in f:
            if post_start:
                line = line.strip()
                pdb_resnum, chain, pikaa, mut_res = line.split()
                mutations.append( [pdb_resnum, chain, pikaa, mut_res] )
            elif line.startswith('start'):
                post_start = True
    return mutations

def find_neighbors(filenum_dir, pdb_path, neighbor_distance = 8.0):
    mutations = read_mutations_resfile(filenum_dir)

    open_filename = pdb_path
    parser = PDBParser(PERMISSIVE=1)

    open_strct = parser.get_structure('Open', open_filename)

    # There should only be one model in PDB file
    num_models = 0
    for model in open_strct.get_models():
        num_models += 1
    assert( num_models == 1 )

    chain_list = [chain.get_id() for chain in open_strct[0].get_chains()]
    neighbors = set()
    for mutation in mutations:
        res_id, chain_id, pikaa, mut_aa = mutation
        mut_chain = str(chain_id)
        try:
            mut_pos = int( res_id )
            mut_insertion_code = ' '
        except ValueError:
            mut_pos = int( res_id[:-1] )
            mut_insertion_code = res_id[-1]

        mut_residue = open_strct[0][mut_chain][(' ', mut_pos, mut_insertion_code)]
        for chain in chain_list:
            for residue in [res.get_id() for res in open_strct[0][chain].get_residues()]:
                try:
                    # Kyle note - might be good to do something else for consistency, since not all residues have CB
                    dist = mut_residue['CB'] - open_strct[0][chain][residue]['CB']
                    if dist < neighbor_distance:
                        neighbors.add( (residue, chain) )
                except KeyError:
                    pass

    return neighbors

#Parses dataset .json file and outputs chain to move and input PDB file directory
def json_parser():
    jsonload = open("data/blank_job_dict.json")
    jsonfile = json.load(jsonload)
    key = sorted(jsonfile.keys())[sge_task_id-1]
    chaintomove = jsonfile[key]["%%chainstomove%%"]
    inputdir = jsonfile[key]['input_file_list'][0]
    
    return chaintomove, inputdir

#Finds neighbors within 8A and adds position and Chain information to a pandas dataframe
def neighbors_list(pdb_filepath, pdb_file):
    neighbors = find_neighbors(pdb_filepath, pdb_file, 8)
    pivotlist = ''
    for i in neighbors:
        string_parse = re.sub("[(),']",'', str(i))
        for s in string_parse.split():
            if s.isdigit():
                pivotlist = pivotlist + s
            else:
                pivotlist = pivotlist + s + ','
                
    ####ADD... turn all data in info dataframe into comma delimited <resnum><chain> pairs
    pivotlist = pivotlist[:-1]
    return pivotlist

#Reads resfile and returns mutation position+chain and type
def resfile_stuff(pdb_filepath):
    resfile = read_mutations_resfile(pdb_filepath)
    position = []
    for i in resfile:
        position.append(i[0])
    return position
    
#Prints CMD input with PDBID, associated mutation, and pivot residues
def bash(chaintomove, inputdir, outputdir):
    #Removes PDB file from path, saves in variable filenum_dir
    inputdir_parse = re.sub("/",' ', str(inputdir))
    data, filenum, pdbtemp = inputdir_parse.split()
    filenum_dir = data + "/" + filenum
    PDBID = pdbtemp[:-4]
    
    #Makes a folder for data dumping
    print 'Making directory %s%s...' %(outputdir, filenum)
    os.mkdir(outputdir + filenum)
  
    #Assigns function output to variables for bash input (pivot_residues, target, resfile_relpath)
    target = resfile_stuff(data_dir)
    pivot_residues = neighbors_list(filenum_dir, inputdir)
    resfile_relpath = os.path.relpath(filenum_dir, outputdir+filenum)
    pdb_relpath = os.path.relpath(inputdir, outputdir+filenum)
    
    targetlist = ''
    for i in target:
        targetlist = targetlist + i + ','
    targetlist = targetlist[:-1]
    
    arg = ['/netapp/home/james.lucas/rosetta_src_2016.08.58479_bundle/main/source/bin/rosetta_scripts.linuxgccrelease',
           '-s',
           pdb_relpath,
           '-parser:protocol',
           'DDG_Test.xml',
           '-ignore_unrecognized_res',
           '-out:path:pdb',
           outputdir + filenum,
           '-parser:script_vars',
           'target=%s' %(targetlist),
           'resfile_relpath=%s' %(resfile_relpath),
           'pivot_residues=%s' %(pivot_residues),
           'chain=%s' %(chaintomove),
           '-nstruct 100'
          ] 
    print 'Working on: %s %s' %(filenum, PDBID)
    
    outfile_path = os.path.join(outputdir+filenum, 'rosetta.out')
    rosetta_outfile = open(outfile_path, 'w')
    print 'Running RosettaScript...'
    rosetta_process = subprocess.Popen(arg, stdout=rosetta_outfile, cwd=outputdir+filenum)
    return_code = rosetta_process.wait()
    print 'Task return code:', return_code, '\n'
    rosetta_outfile.close()    

#Define paths
outputdir = 'output/'

#ACTION!!!
chaintomove, inputdir = json_parser()
bash(chaintomove, inputdir, outp utdir)

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
        
error_out = 'SubmitRun_DDG_Zemu_General.py.e' + job_id + '.' + sge_task_id
output_out = 'SubmitRun_DDG_Zemu_General.py.o' + job_id + '.' + sge_task_id

os.shutil.move(error_out , outputdir + filenum + '/' + error_out)
os.shutil.move(output_out , outputdir + filenum + '/' + output_out)
