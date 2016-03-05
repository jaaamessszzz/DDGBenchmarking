from Bio.PDB import *
import sys
import os
import shutil
import pprint
import klab.cluster_template.parse_settings as parse_settings
import time
import getpass
import json
import re
import shutil
import importlib
import imp
from klab.cluster_template.write_run_file import process as write_run_file
from klab.Reporter import Reporter
import cPickle as pickle
import json

def from_pooja():
    #####
    #This script does 4 things. 1. It find neighbors (within an 8 angstrom shell of the mutation); 2. It generates resfiles and adds a list of neighbors to it; 3. It submits the initial minimization; and 4. It converts pdb numbers to rosetta numbers for the list of neighbor residues and saves it to a file called pivot_res.txt (The list written to this file will later be used for backrub)
    #####

    unique_pdb_num=[]
    for elem in pdb_num:
        if elem not in unique_pdb_num:
            unique_pdb_num.append(elem)
    #print unique_pdb_num

    mutate_res=open("mutate.res", "w")
    repack_res=open("repack.res", "w")
    revert_res=open("revert.res", "w")

    mutate_res.write("NATRO\nstart\n%s %s PIKAA %s\n" %(res_id, chain_id, mut_aa))
    repack_res.write("NATRO\nstart\n")
    revert_res.write("NATRO\nstart\n%s %s PIKAA %s\n" %(res_id, chain_id, wt_aa))

    for i in neighbors:
        repack_res.write(str(i)+"\n")
        if i.split()[0] != res_id:
            mutate_res.write(str(i)+"\n")
            revert_res.write(str(i)+"\n")
    mutate_res.close()
    repack_res.close()
    revert_res.close()
    try:
        os.mkdir("wildtype")
        os.mkdir("mutation")
    except:
        pass
    os.system("cp repack.res revert.res mutation/")
    os.system("cp repack.res mutate.res wildtype/")

    mutation_suffix="_"+str(json_mutation[1:])
    os.system("qsub submit_minimize.py %s %s"%(PDBID[0], mutation_suffix))

    sys.path.insert(0, '/netapp/home/psuresh2/DDG')

    from tools.fs.fsio import read_file
    from tools.rosetta.map_pdb_residues import get_pdb_contents_to_pose_residue_map

    pdb_content = read_file('./%s'%sys.argv[1])
    rosetta_scripts_binary = '/netapp/home/thompson/Rosetta/rosetta_2014.22/source/bin/rosetta_scripts.linuxgccrelease'
    success, result = get_pdb_contents_to_pose_residue_map(pdb_content, rosetta_scripts_binary, pdb_id = None, extra_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res')
    #pprint.pprint(result)

    rosetta_num=[]
    for i in range(len(unique_pdb_num)):
        #print result[str(unique_pdb_num[i])]['pose_residue_id']
        rosetta_num.append(result[str(unique_pdb_num[i])]['pose_residue_id'])

    #print rosetta_num

    pivot_residues=open("pivot_res.txt", "w")
    res=str(" ".join([str(x) for x in rosetta_num]))
    pivot_residues.write(res)

    minimize_pdb=open("pdblist.txt", "w")
    minimize_pdb.write("%s_0001.pdb\n%s%s_0001.pdb"%(str(PDBID[0]), str(PDBID[0]), str(mutation_suffix)))

def find_neighbors(data_dir, pdb_path, neighbor_distance = 8.0):
    resfile = os.path.join(data_dir, 'mutations.resfile')
    mutations = []
    with open(resfile, 'r') as f:
        post_start = False
        for line in f:
            if post_start:
                line = line.strip()
                pdb_resnum, chain, pikaa, mut_res = line.split()
                mutations.append( [pdb_resnum, chain, pikaa, mut_res] )
            elif line.startswith('start'):
                post_start = True

    open_filename = pdb_path
    parser = PDBParser(PERMISSIVE=1)

    open_strct = parser.get_structure('Open', open_filename)

    # There should only be one model in PDB file
    num_models = 0
    for model in open_strct.get_models():
        num_models += 1
    assert( num_models == 1 )

    chain_list = [chain.get_id() for chain in open_strct[0].get_chains()]
    neighbors = set()
    for mutation in mutations:
        res_id, chain_id, pikaa, mut_aa = mutation
        ### print res_id, chain_id, pikaa, mut_aa
        mut_chain = str(chain_id)
        try:
            mut_pos = int( res_id )
            mut_insertion_code = ' '
        except ValueError:
            mut_pos = int( res_id[:-1] )
            mut_insertion_code = res_id[-1]

        mut_residue = open_strct[0][mut_chain][(' ', mut_pos, mut_insertion_code)]
        for chain in chain_list:
            for residue in [res.get_id() for res in open_strct[0][chain].get_residues()]:
                try:
                    # Kyle note - might be good to do something else for consistency, since not all residues have CB
                    dist = mut_residue['CB'] - open_strct[0][chain][residue]['CB']
                    if dist < neighbor_distance:
                        neighbors.add( (residue, chain) )
                except KeyError:
                    pass

    return neighbors

def write_repack_resfile(resfile_path, neighbors):
    with open(resfile_path, 'w') as f:
        f.write('NATRO\nstart\n')
        for neighbor in sorted(list(neighbors)):
            res_tup, chain = neighbor
            het_flag, res_num, insertion_code = res_tup
            insertion_code = insertion_code.strip()
            f.write( '%d%s %s NATAA\n' % (res_num, insertion_code, chain) )

def neighbors_to_rosetta_num(neigbors, data_dir):
    numbering_map_path = os.path.join(data_dir, 'pdb2rosetta.resmap.json')
    with open(numbering_map_path, 'r') as f:
        numbering_map = json.load(f)
    numbering_map = {
        (key.strip().split()[0], key.strip().split()[1]) : value
        for key, value in numbering_map.iteritems()
    }
    return_nums = set()
    for tup, chain in neighbors:
        het_flag, pdb_num, insertion_code = tup
        return_nums.add( numbering_map[ (chain, '%d%s' % (pdb_num, insertion_code.strip()) ) ] )
    return sorted( list(return_nums) )

if __name__ == '__main__':
    assert( len(sys.argv) > 1 )
    cfg = imp.load_source('cfg', sys.argv[1])
    assert( os.path.isdir(sys.argv[2]) ) # Output directory

    prediction_set_id = cfg.prediction_set_id

    output_dir = sys.argv[2]
    data_dir_path = os.path.join(output_dir, 'data')
    job_dict_path = os.path.join(data_dir_path, 'blank_job_dict.pickle')
    settings_path = os.path.join(data_dir_path, 'settings.pickle')

    with open(job_dict_path, 'r') as f:
        job_dict = pickle.load(f)
    with open(settings_path, 'r') as f:
        settings = pickle.load(f)

    r = Reporter('processing prediction ids', entries = 'prediction ids')
    r.set_total_count( len(job_dict) )
    for prediction_id in job_dict:
        assert( len(job_dict[prediction_id]['input_file_list']) == 1 )
        pdb = job_dict[prediction_id]['input_file_list'][0]
        pdb_path = os.path.join(output_dir, pdb)
        data_dir = os.path.join( os.path.join(output_dir, 'data'), '%d' % prediction_id )
        neighbors = find_neighbors(data_dir, pdb_path)
        repack_resfile_path = os.path.join(data_dir, 'repack.resfile')
        repack_resfile_relpath = os.path.relpath(repack_resfile_path, output_dir)
        pivot_residues = ' '.join( ['%d' % d for d in neighbors_to_rosetta_num(neighbors, data_dir)] ).strip()
        job_dict[prediction_id]['-backrub:pivot_residues'] = pivot_residues
        print job_dict[prediction_id]
        write_repack_resfile(repack_resfile_path, neighbors)
        break
        r.increment_report()
    r.done()

    settings['output_dir'] = output_dir
    settings['tasks_per_process'] = 1
    settings['mem_free'] = '1.6G'
    settings['rosetta_args_list'].extend( [
        '-nstruct 50',
        '-ignore_unrecognized_res',
        '-ignore_zero_occupancy false',
        '-ex1', '-ex2',
        '-extrachi_cutoff 0',
        '-out:prefix bkrb_',
        '-mute core.io.pdb.file_data',
        '-backrub:ntrials 50000',
        '-mc_kt %.1f' % cfg.backrub_temp,
    ] )
    settings['appname'] = 'backrub'

    # print job_dict
    # print settings
    write_run_file(settings, database_run = False, job_dict = job_dict)
    print 'Job files written to directory:', os.path.abspath(output_dir)
