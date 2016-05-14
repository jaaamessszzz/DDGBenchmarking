#!/bin/bash

NUMPROC=$(nproc)

# Function to wait until less jobs are running than number of processors
wait_for_free_processor() {
    joblist=($(jobs -p))
    # echo "Running jobs: ${#joblist[*]}"
    while (( ${#joblist[*]} >= $NUMPROC ))
    do
    # echo "${#joblist[*]} jobs running"
        sleep 1
        joblist=($(jobs -p))
    done
}

cd /kortemmelab/home/james.lucas/160412-kyleb_jl-brub-rscr-v2/DDG_Zemu_v2_output-Sum_DDG_only

#Analyze PDBs from one directory
nice python ../../DDGBenchmarks_Test/UtilityScripts/RMSD_X_Calculator.py 68417 &

#Iterate over subdirectories within whatever starting directory
#for i in * ; do
#    if [ -d "$i" ]; then
#	nice python ../../DDGBenchmarks_Test/UtilityScripts/RMSD_X_Calculator.py $i &
#	wait_for_free_processor
#    fi
#done
#wait
