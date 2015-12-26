#!/bin/bash
for i in {1..15}; do
    echo $i
    date
    python analyze_old_ubq_runs.py /dbscratch/kyleb/tmp/cluster_run/151115-kyleb_rescore_ddg_monomer-$i
    date
    sleep 10
done




