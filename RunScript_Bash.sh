#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_rt=24:00:00
##$ -t 1-200
#$ -l arch=linux-x64
#$ -l mem_free=1.5G

~/rosetta_src_2016.08.58479_bundle/main/source/bin/rosetta_scripts.linuxgccrelease -s ~/107L.pdb -parser:protocol ~/DDG_Test.xml -ignore_unrecognized_res -out:path:pdb ~/Test_Output -nstruct 100