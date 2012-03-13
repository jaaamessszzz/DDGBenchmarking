numprocs=30
jobname=$1

#exe & database
exepath="/ifs1/home/stranges/rosetta/mini_optE_fixes/build/src/release/linux/2.6/64/x86/gcc/mpi/surface_optE_parallel.linuxgccrelease"
dbpath="/ifs1/home/stranges/rosetta/minirosetta_database/"

#for sequence recovery
#pdblist="/ifs1/home/stranges/scr/interface_optE_inputs/nataa_recovery_pdbs.list"
#raslist="/ifs1/home/stranges/scr/interface_optE_inputs/ras_xtals.list"
ddmilist="/ifs1/home/stranges/scr/interface_optE_inputs/ddmi_xtals_list.txt"

#what mut data to use?
#ddGbinddata="/ifs1/scr/stranges/interface_optE_inputs/ras_lists.ddg.input"
#ddGbinddata="/ifs1/scr/stranges/interface_optE_inputs/bder_lists.ddg.input"
#ddGbinddata="/ifs1/scr/stranges/interface_optE_inputs/bder_no1dan_lists.ddg.input"
ddGbinddata="/ifs1/scr/stranges/interface_optE_inputs/bder_best_lists.ddg.input"

#tagfile for task making
#tagfile="/ifs1/scr/stranges/interface_optE/017c.pNATAA_seqpropensity_sametask/prot_interface_design_task.xml"
tagfile="/ifs1/scr/stranges/interface_optE/017d.pNATAA_seqpropensity_sametask/prot_interface_design_task.xml"
const_tagfile="/ifs1/scr/stranges/interface_optE/017d.pNATAA_seqpropensity_sametask/const_task.xml"

logfile="$jobname.run.log"

bsubcommand="bsub -q day -n $numprocs -J $jobname -o $logfile -a mvapich mpirun"

cmd="$bsubcommand $exepath -database $dbpath \
-s $ddmilist \
-optE:optimize_nat_aa true \
-optE:optimize_ddG_bind_correlation $ddGbinddata \
-optE:n_design_cycles 15 \
-optE:mpi_weight_minimization true \
-optE:fit_reference_energies_to_aa_profile_recovery true \
-optE:repeat_swarm_optimization_until_fitness_improves true \
-optE:optimize_starting_free_weights true -optE:number_of_swarm_particles 25 -optE:number_of_swarm_cycles 10 \
-optE::component_weights fitness_function_component_weights.txt \
-optE:free free_wts.txt -optE:fixed fixed_wts.txt \
-optE:parse_tagfile $tagfile \
-optE:constant_logic_taskops_file $const_tagfile \
-ignore_unrecognized_res -no_optH false -skip_set_reasonable_fold_tree -no_his_his_pairE \
-ex1 -ex2 -linmem_ig 5 -extrachi_cutoff 1 -options:user \
-mute core.io core.pack core.scoring core.conformation ig.node ig.bgnode ig.edge ig.bgedge protocols.optimize_weights.OptEMultifunc core.optimization.LineMinimizer \
-exclude_badrep_ddGs 5.0"

# clean up old directories, make fresh ones
rm -rf workdir_* weightdir logdir

i=0
while [ $i -lt $numprocs ]
do
  mkdir "workdir_$i"
  i=$(($i+1))
done
mkdir weightdir
mkdir logdir

# print and execute command
echo $cmd
$cmd

#-mute core.io core.pack core.scoring core.conformation ig.node ig.bgnode ig.edge ig.bgedge protocols.optimize_weights.OptEMultifunc core.optimization.LineMinimizer \

#-ex1 -ex2 -linmem_ig 5 -extrachi_cutoff 1 -options:user \
#-optE:constant_logic_taskops_file $const_tagfile \
