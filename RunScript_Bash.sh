#!/usr/bin/bash
#$ -S /usr/bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_rt=24:00:00
##$ -t 1-200
#$ -l arch=linux-x64
#$ -l mem_free=1.5G

./rosetta_scripts.macosclangrelease -s ~/107L.pdb -parser:protocol ~/DDG_Test.xml -ignore_unrecognized_res -out:path:pdb ~/Test_Output -nstruct 100