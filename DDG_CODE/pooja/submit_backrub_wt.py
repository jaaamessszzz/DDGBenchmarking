#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=48:00:00
#$ -t 1-100
#$ -l arch=linux-x64
#$ -l mem_free=4.0G
#$ -l netapp=1G

import socket
import sys

# Print the python/host information as soon as possible
print "Python:", sys.version
print "Host:", socket.gethostname()

import datetime
import os
import subprocess
import time

iteration = 1
if os.environ.has_key("SGE_TASK_ID"):
        iteration = int(os.environ["SGE_TASK_ID"])

PDBFILE=str(sys.argv[1])

piv_res_file=open("pivot_res.txt", "r")
pivot_residues=piv_res_file.readline()

backrub_command=str("/netapp/home/thompson/Rosetta/rosetta_2014.22/source/bin/backrub.linuxgccrelease -database /netapp/home/thompson/Rosetta/rosetta_2014.22/database -s %s -ignore_unrecognized_res -ignore_zero_occupancy false -ex1 -ex2 -pivot_residues %s -extrachi_cutoff 0 -out:prefix bkrb_%04i_ -mute core.io.pdb.file_data -backrub:ntrials 50000 -mc_kt 1.2 -overwrite") %(PDBFILE, pivot_residues, iteration)

print backrub_command

subprocess.call(backrub_command, shellx=True)
###############repack and mutate
from glob import glob
file_list=glob("*last.pdb")

iteration = 1
if os.environ.has_key("SGE_TASK_ID"):
        iteration = int(os.environ["SGE_TASK_ID"])

repack_command=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile repack.res -out:suffix _repack -ex1 -ex2 -multi_cool_annealer 6 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite") %file_list[iteration-1]

mutate_command=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile mutate.res -out:suffix _mutate -ex1 -ex2 -multi_cool_annealer 6 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite") %file_list[iteration-1]

print repack_command
print mutate_command

subprocess.call(repack_command, shell=True)
subprocess.call(mutate_command, shell=True)
