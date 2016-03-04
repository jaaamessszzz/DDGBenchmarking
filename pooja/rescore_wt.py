#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=48:00:00
#$ -t 1-100
#$ -l arch=linux-x64
#$ -l mem_free=4.0
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

from glob import glob
file_list_mut_com=glob("*_0001_0001_0001_last_mutate_0001_0001.pdb")
file_list_rep_com=glob("*_0001_0001_0001_last_repack_0001_0001.pdb")
file_list_mut_mon=glob("*_0001_0001_0001_last_mutate_0001_0001_monomer.pdb")
file_list_rep_mon=glob("*_0001_0001_0001_last_repack_0001_0001_monomer.pdb")
file_list_mut_par=glob("*_0001_0001_0001_last_mutate_0001_0001_partner.pdb")
file_list_rep_par=glob("*_0001_0001_0001_last_repack_0001_0001_partner.pdb")

iteration = 1
if os.environ.has_key("SGE_TASK_ID"):
        iteration = int(os.environ["SGE_TASK_ID"])

rescore_mut_com=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile natro.res -out:suffix _natro -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite")%file_list_mut_com[iteration-1]

rescore_rep_com=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile natro.res -out:suffix _natro -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite")%file_list_rep_com[iteration-1]

rescore_mut_mon=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile natro.res -out:suffix _natro -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite")%file_list_mut_mon[iteration-1]

rescore_rep_mon=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile natro.res -out:suffix _natro -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite")%file_list_rep_mon[iteration-1]

rescore_mut_par=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile natro.res -out:suffix _natro -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite")%file_list_mut_par[iteration-1]

rescore_rep_par=str("/netapp/home/thompson/Rosetta/rosetta_2014.35/source/bin/fixbb.default.linuxgccrelease -s %s -database /netapp/home/thompson/Rosetta/rosetta_2014.35/database -resfile natro.res -out:suffix _natro -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true -overwrite")%file_list_rep_par[iteration-1]


print rescore_mut_com
print rescore_rep_com
print rescore_mut_mon
print rescore_rep_mon
print rescore_mut_par
print rescore_rep_par

subprocess.call(rescore_mut_com, shell=True)
subprocess.call(rescore_rep_com, shell=True)
subprocess.call(rescore_mut_mon, shell=True)
subprocess.call(rescore_rep_mon, shell=True)
subprocess.call(rescore_mut_par, shell=True)
subprocess.call(rescore_rep_par, shell=True)
