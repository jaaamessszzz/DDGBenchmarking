#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-707
#$ -l arch=linux-x64
#$ -l mem_free=1.1G
#$ -l netapp=1G,scratch=1G

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

def read_mutations_resfile(data_dir):
    resfile = os.path.join(data_dir, 'mutations.resfile')
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

def find_neighbors(data_dir, pdb_path, neighbor_distance = 8.0):
    mutations = read_mutations_resfile(data_dir)

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
    
####asdf = open("~/ddg-zemu-runs/160315-kyleb_james-backrub-rscript/data/blank_job_dict.json")
    asdf = open("/netapp/home/james.lucas/160315-kyleb_james-backrub-rscript/data/blank_job_dict.json")
    jsonfile = json.load(asdf)
    key = sorted(jsonfile.keys())[sge_task_id-1]
    chaintomove = jsonfile[key]["%%chainstomove%%"]
    directory = jsonfile[key]['input_file_list'][0]
    
    return chaintomove, directory

#Finds neighbors within 8A and adds position and Chain information to a pandas dataframe
def neighbors_list(pdb_filepath, pdb_file):
    neighbors = find_neighbors("/netapp/home/james.lucas/160315-kyleb_james-backrub-rscript/"+pdb_filepath, "/netapp/home/james.lucas/160315-kyleb_james-backrub-rscript/"+pdb_file, 8)
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
    resfile = read_mutations_resfile("/netapp/home/james.lucas/160315-kyleb_james-backrub-rscript/"+pdb_filepath)
    for i in resfile:
        position = i[0]
        chain = i[1]
        mut_to = i[3]
    return position, chain, mut_to
    
#Prints CMD input with PDBID, associated mutation, and pivot residues
def bash(chaintomove, pdb_file):
    #Working Directory
    workingdir = '/netapp/home/james.lucas/DDG_Zemu_Backrub_Output/'
    #Removes PDB file from path, saves in variable data_dir
    pdb_file_parse = re.sub("/",' ', str(pdb_file))
    data, filenum, pdbtemp = pdb_file_parse.split()
    data_dir = data + "/" + filenum
    PDBID = pdbtemp[:-4]
    
    #Makes a folder for data dumping
    print 'Making directory %s%s...' %(workingdir, filenum)
    os.mkdir(workingdir + filenum)
        
    #Dictionary: 1- to 3-letter code
    res_dict = {
        'A':'ALA',
        'C':'CYS',
        'D':'ASP',
        'E':'GLU',
        'F':'PHE',
        'G':'GLY',
        'H':'HIS',
        'I':'ILE',
        'K':'LYS',
        'L':'LEU',
        'M':'MET',
        'N':'ASN',
        'P':'PRO',
        'Q':'GLN',
        'R':'ARG',
        'S':'SER',
        'T':'THR',
        'V':'VAL',
        'W':'TRP',
        'Y':'TYR'
    }
        
    #Assigns function output to variables for bash input (pivot_residues, target, new_res)
    target, chain, new_res_one = resfile_stuff(data_dir)
    new_res_three = res_dict[new_res_one]
    pivot_residues = neighbors_list(data_dir, pdb_file)
        
####All the variables and stuff for printing out the bash script
    
#Cluster Usage - Minimization
    arg = ['/netapp/home/james.lucas/rosetta_src_2016.08.58479_bundle/main/source/bin/rosetta_scripts.linuxgccrelease',
           '-s',
           '/netapp/home/james.lucas/160315-kyleb_james-backrub-rscript/%s/%s.pdb' %(data_dir, PDBID),
           '-parser:protocol',
           '/netapp/home/james.lucas/DDG_Test.xml',
           '-ignore_unrecognized_res',
           '-out:path:pdb',
           workingdir + filenum,
           '-parser:script_vars',
           'target=%s' %(target),
           'new_res=%s' %(new_res_three),
           'pivot_residues=%s' %(pivot_residues),
           'chain=%s' %(chain),
           '-nstruct 100'
          ] 
    print 'Working on: %s %s%s' %(PDBID, target, new_res_three)
    
    outfile_path = os.path.join(workingdir+filenum, 'rosetta.out')
    rosetta_outfile = open(outfile_path, 'w')
    print 'Running RosettaScript...'
    rosetta_process = subprocess.Popen(arg, stdout=rosetta_outfile, cwd=workingdir+filenum)
    return_code = rosetta_process.wait()
    print 'Task return code:', return_code, '\n'
    rosetta_outfile.close()
    
    #subprocess.call(arg)    
    
#ACTION!!!
chaintomove, pdb_file = json_parser()
bash(chaintomove, pdb_file)

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
