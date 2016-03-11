#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_rt=24:00:00
#$ -t 1
#$ -l arch=linux-x64
#$ -l mem_free=1.5G

import socket
import sys
import os
import subprocess
import shutil

print "Python:", sys.version
print "Host:", socket.gethostname()

sge_task_id = 1
if os.environ.has_key("SGE_TASK_ID"):
	sge_task_id = os.environ["SGE_TASK_ID"]
id = int(sge_task_id)

cwd = os.getcwd()

rosettasripts_path = "/netapp/home/james.lucas/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease"
xml = "/netapp/home/james.lucas/DDG_Test.xml"
out_path = "~/test_results"
database_path = "/netapp/home/james.lucas/"

#./rosetta_scripts.macosclangrelease -s ~/Kortemme_Rotation/107L.pdb -parser:protocol ~/Kortemme_Rotation/DDG_Test.xml -overwrite -ignore_unrecognized_res -out:path:pdb ~/Kortemme_Rotation/Output -nstruct 100

backrub_args = [
	backrub_path,
	"-database %s" % database_path,
    xml,
	"-s %s" % pdb_file,
	"-nstruct 100",
	"-ignore_unrecognized_res",
    "-out:path:pdb %s" % outpath
	"-out:prefix %s_" %id,
]

print backrub_args

try:
	os.mkdir('test_results')
except:
	pass

new_id = '000' + str(id)
new_id = new_id[-4:]
pdb_file_base = pdb_file.split('/')[-1].split('.pdb')[0]
new_file = pdb_file_base+ '_' +new_id + '_last.pdb' 
old_file = str(id) + '_' + pdb_file_base + '_0001_last.pdb'

os.chdir('test_results')

if os.path.isfile(new_file) == False:
    outfile = open(name+'.log', 'w')
    process = subprocess.Popen(backrub_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True)
    returncode = process.wait()
    shutil.move(old_file, new_file)
    outfile.close()

elif os.path.isfile(new_file) == True:
	print 'EXIT: file exists'
