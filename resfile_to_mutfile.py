#!/usr/bin/python

## simple script to convert a standard Rosetta resfile into a mutfile for use with Liz' ddG protocol
## input: a resfile and the PDB it refers to
## output: an individual mutfile for each possible mutation indicated in the resfile, to enable
##         parallel computation of the results, and a script to run all of them
##
## NOTE: - for now, this is ignoring EX [ARO] flags given in the resfile
##       - all mutations are considered independent, this script will not generate mutfiles for multi-site mutations
##
## author: Amelie Stein (amelie.stein@ucsf.edu)

import operator
import os
from rosetta import *


##
## parameters you may want to change
##
run_cmd = "python" ## replace by qsub for running the individual jobs in parallel on the cluster
ddG_protocol = "ddg_protocol13_constraints_mutfile.py" ## the protocol wrapper
## if we need parameter files etc., include them here


## function to write a mutfile, translated from each possibility in the resfile
def write_mutfile(wt_res, pos, target_res, filename):
    file = open(filename, "w") ## should include the subdirectory (if relevant) 
    file.write('total 1\n1\n')
    file.write('%s %s %s\n' % (wt_res, pos, target_res))
    file.close()
#-



def main():
        if (len(sys.argv) != 3):
                print "Error -- you need to provide the resfile and the [ideally renumbered] PDB file\n";
                return
        #-
        resfile = sys.argv[1]
        pdb_file = sys.argv[2] # can be fully qualified path

## generate a script to run all the mutations in this resfile
        run_file_name = "run_mutations_from_"+resfile+".sh"
        run_file = open(run_file_name, "w")

## also generate a mapping between the chain+seqres and the Rosetta position?

## later: also generate a script that runs all the evaluation commands -- TODO
        

## we'd need to check for .params files here if a ligand is present -- CURRENTLY NOT IMPLEMENTED
        rosetta.init()
        p = Pose()
        pose_from_pdb(p, pdb_file)

        d = pdb_file.split('/')
        pdb_id = d.pop()
        #pdb_path = "/".join(d) ## needed?
        if (pdb_id.endswith(".pdb")):
            pdb_id = pdb_id[0:len(pdb_id)-4]
        #-

## write the constraints file -- only needed once for each PDB
        print "[INFO] writing the constraint file for "+pdb_id+"..."
        constraints_file_name = pdb_id+".constraints"
        constraints_file = open(constraints_file_name, "w")
        for i in range(1,p.total_residue()+1):
            res_i = p.residue(i)
            if res_i.is_protein(): ## ions and ligands don't necessarily have CAs and will cause a crash
                for j in range(i+1,p.total_residue()+1):
                    res_j = p.residue(j)
                    if res_j.is_protein():
                        dist = res_i.xyz("CA").distance(res_j.xyz("CA"))
                        if dist < 9.0:
                            constraints_file.write("AtomPair CA "+str(i)+" CA "+str(j)+" HARMONIC "+str(dist)+" 0.5\n")
        #------
        constraints_file.close()


## base directory for all the mutations and files -- subdirectory of the current location
        os.popen('mkdir ' + pdb_id) 

## now create a mutfile for each mutation        
        rf = open(resfile, 'r')
        for line in rf:
            ## currently we can only handle PIKAA, others should probably be implemented eventually -- cf. http://www.rosettacommons.org/manuals/rosetta3_user_guide/file_resfiles.html
            #if line.find("PIKAA") != None:
            if "PIKAA" in line:
                #print "handling "+line ## debugging
                [pos, chain, cmd, target_res_set] = line.split()
                
            ## lookup matching Rosetta numbering and create a subdirectory for the position
                rosetta_pos = p.pdb_info().pdb2pose(chain, int(pos))
                print "[INFO] Position "+pos+" in chain "+chain+" corresponds to "+str(rosetta_pos)+" in Rosetta numbering"
                pos_dir = pdb_id+"/"+str(rosetta_pos)
                os.popen('mkdir ' + pos_dir)
                run_file.write('cd %s\n' % (pos_dir))

                ## copy/link the structure, the constraint file and the script into this directory
                if pdb_file.startswith("/"):
                    os.popen('ln -s '+pdb_file+' %s' % (pos_dir))
                else:
                    os.popen('ln -s ../../'+pdb_file+' %s' % (pos_dir))
                #-
                os.popen('ln -s ../../'+constraints_file_name+' %s' % (pos_dir))
                os.popen('ln -s ../../'+ddG_protocol+' %s' % (pos_dir))
                
                wt_res = p.residue(rosetta_pos).name1() ## mutfiles need to know the wt residue
                if wt_res not in target_res_set:
                    target_res_set += wt_res ## we also need to generate structures for the WT to calculate the ddG
                for target_res in target_res_set:
                    mutfile = wt_res+str(rosetta_pos)+target_res+".mut"
                    # print "trying to generate "+mutfile ## debugging
                    write_mutfile(wt_res, rosetta_pos, target_res, pos_dir+"/"+mutfile)
                    run_file.write('%s %s %s %s %s %s\n' % (run_cmd, ddG_protocol, pdb_id, rosetta_pos, constraints_file_name, mutfile))
                #-
                run_file.write('cd -\n')

        #--
        
        run_file.close()

if __name__ == "__main__":
        main()
