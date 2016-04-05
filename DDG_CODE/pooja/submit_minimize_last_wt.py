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

pdblist=open("pdblist.txt", "w")
from glob import glob
file_list_mut=glob("*last_mutate_0001.pdb")
file_list_rep=glob("*last_repack_0001.pdb")
file_list_mut.sort()
file_list_rep.sort()
for i in range(len(file_list_mut)):
    pdblist.write(str(file_list_mut[i])+"\n")
    pdblist.write(str(file_list_rep[i])+"\n")

minimize_command=str("/netapp/home/klabqb3backrub/r57200/minimize_with_cst.static.linuxgccrelease -database /netapp/home/klabqb3backrub/r57200/rosetta_database -in:file:l pdblist.txt -in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -ddg::harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst -ddg::sc_min_only false -ignore_zero_occupancy false -overwrite")

subprocess.call(minimize_command, shell=True)
