#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=48:00:00
#$ -t 1
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

PDBID=str(sys.argv[1])
mutation=str(sys.argv[2])

mut_fixbb_command=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s.pdb -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile mutate.res -out:suffix %s -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite") %(PDBID, mutation)

repack_fixbb_command=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s.pdb -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile repack.res -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite") %PDBID

print mut_fixbb_command
subprocess.call(mut_fixbb_command, shell=True)

print repack_fixbb_command
subprocess.call(repack_fixbb_command, shell=True)

minimize_command=str("/netapp/home/klabqb3backrub/r57200/minimize_with_cst.static.linuxgccrelease -database /netapp/home/klabqb3backrub/r57200/rosetta_database -in:file:l pdblist.txt -in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -ddg::harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst -ddg::sc_min_only false -ignore_zero_occupancy false -overwrite")

print minimize_command
subprocess.call(minimize_command, shell=True)
