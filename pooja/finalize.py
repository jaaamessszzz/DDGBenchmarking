from Bio.PDB import *
import sys
import os
import shutil
import pprint
import klab.cluster_template.cluster_template as cluster_template
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
import copy

def read_mutations_resfile(data_dir):
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
    return mutations

def find_neighbors(data_dir, pdb_path, neighbor_distance = 8.0):
    mutations = read_mutations_resfile(data_dir)

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

def write_repack_resfile(resfile_path, neighbors, data_dir = None):
    if data_dir:
        # Parse mutation resfile and add to result resfile
        mutations = {(res_id, chain_id) : mut_aa for res_id, chain_id, pikaa, mut_aa in read_mutations_resfile(data_dir)}
    else:
        mutations = {}

    with open(resfile_path, 'w') as f:
        f.write('NATRO\nstart\n')
        for neighbor in sorted(list(neighbors)):
            res_tup, chain = neighbor
            het_flag, res_num, insertion_code = res_tup
            insertion_code = insertion_code.strip()
            res_id = '%d%s' % (res_num, insertion_code)
            if (res_id, chain) in mutations:
                f.write( '%s %s PIKAA %s\n' % (res_id, chain, mutations[(res_id, chain)]) )
            else:
                f.write( '%s %s NATAA\n' % (res_id, chain) )
    if len(mutations) > 1:
        print resfile_path, mutations
        sys.exit(1)

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

    # This uses the version of Rosetta from your cluster template settings file
    ppi_api_lib = imp.load_source('ppi_api', '../ddglib/ppi_api.py')
    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']
    ppi_api = ppi_api_lib.get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database')

    with open(job_dict_path, 'r') as f:
        original_job_dict = pickle.load(f)
    with open(settings_path, 'r') as f:
        settings = pickle.load(f)
    settings['output_dir'] = output_dir
    settings['tasks_per_process'] = 1
    settings['mem_free'] = '1.6G'
    num_steps = 6 # Preminimization, backrub, repack wt, mutate, minimize repack wt, minimize mutant
    settings = cluster_template.convert_list_arguments_to_list(settings, num_steps)
    settings['appname_list'] = []
    inner_jobs = 2 ### TMP, should be 100
    settings['numjobs'] = settings['numjobs'] * inner_jobs

    job_dicts = []

    generic_min_args = [
        '-in:file:fullatom',
        '-ignore_unrecognized_res',
        '-fa_max_dis 9.0',
        '-ddg::harmonic_ca_tether 0.5',
        '-ddg::constraint_weight 1.0',
        '-ddg::out_pdb_prefix min_cst',
        '-ddg::sc_min_only false',
        '-ignore_zero_occupancy false',
    ]

    # Preminimization step
    print 'Adding minimization step\n'
    job_dict = {}
    settings['rosetta_args_list_list'][len(job_dicts)].extend(
        generic_min_args
    )
    for prediction_id in original_job_dict:
        for n_struct in xrange(1, inner_jobs + 1):
            inner_prediction_id = '%s/%d-%d-minimize' % (str(prediction_id), n_struct, len(job_dicts))
            job_dict[inner_prediction_id] = copy.deepcopy( original_job_dict[prediction_id] )
        break ### TMP
    settings['appname_list'].append( 'minimize_with_cst' )
    job_dicts.append( job_dict )

    # Backrub step
    job_dict = {}
    repack_resfile_relpaths = {} # For later step
    repack_mutate_resfile_relpaths = {} # For later step
    pdb_ids_from_pred_ids = {} # For later step
    r = Reporter('processing prediction ids for backrub information', entries = 'prediction ids')
    r.set_total_count( len(original_job_dict) )
    for prediction_id in original_job_dict:
        pdb = job_dicts[0]['%s/%d-%d-minimize' % (str(prediction_id), 1, 0)]['input_file_list'][0]
        pdb_id = os.path.basename(pdb).split('.')[0].strip()
        pdb_ids_from_pred_ids[prediction_id] = pdb_id
        pdb_path = os.path.join(output_dir, pdb)
        data_dir = os.path.join( os.path.join(output_dir, 'data'), '%d' % prediction_id )
        neighbors = find_neighbors(data_dir, pdb_path)
        repack_resfile_path = os.path.join(data_dir, 'repack.resfile')
        repack_resfile_relpath = os.path.relpath(repack_resfile_path, output_dir)
        repack_resfile_relpaths[prediction_id] = repack_resfile_relpath
        pivot_residues = ' '.join( ['%d' % d for d in neighbors_to_rosetta_num(neighbors, data_dir)] ).strip()
        write_repack_resfile(repack_resfile_path, neighbors)
        repack_mutate_resfile_path = os.path.join(data_dir, 'repack_and_mutate.resfile')
        repack_mutate_resfile_relpath = os.path.relpath(repack_mutate_resfile_path, output_dir)
        repack_mutate_resfile_relpaths[prediction_id] = repack_mutate_resfile_relpath
        write_repack_resfile(repack_mutate_resfile_path, neighbors, data_dir = data_dir)
        for n_struct in xrange(1, inner_jobs + 1):
            inner_prediction_id = '%s/%d-%d-backrub' % (str(prediction_id), n_struct, len(job_dicts))
            minimize_inner_prediction_id = '%s/%d-%d-minimize' % (str(prediction_id), n_struct, len(job_dicts)-1)

            inner_prediction_min_dir = os.path.join(output_dir, minimize_inner_prediction_id)
            pdb_id = os.path.basename(pdb).split('.')[0].strip()
            minimized_pdb = os.path.relpath(os.path.join(inner_prediction_min_dir, 'min_cst.%s_0001.pdb.gz' % pdb_id), output_dir)
            job_dict[inner_prediction_id] = {}
            job_dict[inner_prediction_id]['input_file_list'] = [
                minimized_pdb
            ]

            job_dict[inner_prediction_id]['-backrub:pivot_residues'] = pivot_residues
        r.increment_report()
        break ### TMP
    r.done()

    settings['rosetta_args_list_list'][len(job_dicts)].extend( [
        '-nstruct 1',
        '-ignore_unrecognized_res',
        '-ignore_zero_occupancy false',
        '-ex1', '-ex2',
        '-extrachi_cutoff 0',
        '-out:prefix bkrb_',
        '-mute core.io.pdb.file_data',
        '-backrub:ntrials 1', ### TMP for testing - should really be 50000!
        '-mc_kt %.1f' % cfg.backrub_temp,
    ] )
    settings['appname_list'].append( 'backrub' )
    job_dicts.append( job_dict )

    repack_generic_args = [
        '-ex1', '-ex2',
        '-multi_cool_annealer 6',
        '-ignore_zero_occupancy false',
        '-ignore_unrecognized_res true',
    ]

    # Repack wildtype step
    job_dict = {}
    for prediction_id in original_job_dict:
        data_dir = os.path.join( os.path.join(output_dir, 'data'), '%d' % prediction_id )
        pdb_id = pdb_ids_from_pred_ids[prediction_id]
        for n_struct in xrange(1, inner_jobs + 1):
            inner_prediction_id = '%s/%d-%d-repack_wt' % (str(prediction_id), n_struct, len(job_dicts))
            backrub_inner_prediction_id = '%s/%d-%d-backrub' % (str(prediction_id), n_struct, len(job_dicts)-1)

            inner_prediction_backrub_dir = os.path.join(output_dir, backrub_inner_prediction_id)
            backrubed_pdb = os.path.relpath(os.path.join(inner_prediction_backrub_dir, 'bkrb_min_cst.%s_0001_0001_last.pdb.gz' % pdb_id), output_dir)
            job_dict[inner_prediction_id] = {}
            job_dict[inner_prediction_id]['input_file_list'] = [
                backrubed_pdb
            ]
            job_dict[inner_prediction_id]['-resfile'] = repack_resfile_relpaths[prediction_id]
        break ### TMP

    settings['rosetta_args_list_list'][len(job_dicts)].extend(
        repack_generic_args
    )
    settings['rosetta_args_list_list'][len(job_dicts)].extend( [
        '-out:prefix repack-wt_',
    ] )
    settings['appname_list'].append( 'fixbb' )
    job_dicts.append( job_dict )

    # Mutant repack step
    job_dict = {}
    for prediction_id in original_job_dict:
        data_dir = os.path.join( os.path.join(output_dir, 'data'), '%d' % prediction_id )
        pdb_id = pdb_ids_from_pred_ids[prediction_id]
        for n_struct in xrange(1, inner_jobs + 1):
            inner_prediction_id = '%s/%d-%d-mutate' % (str(prediction_id), n_struct, len(job_dicts))
            backrub_inner_prediction_id = '%s/%d-%d-backrub' % (str(prediction_id), n_struct, len(job_dicts)-2)

            inner_prediction_backrub_dir = os.path.join(output_dir, backrub_inner_prediction_id)
            backrubed_pdb = os.path.relpath(os.path.join(inner_prediction_backrub_dir, 'bkrb_min_cst.%s_0001_0001_last.pdb.gz' % pdb_id), output_dir)
            job_dict[inner_prediction_id] = {}
            job_dict[inner_prediction_id]['input_file_list'] = [
                backrubed_pdb
            ]
            job_dict[inner_prediction_id]['-resfile'] = repack_mutate_resfile_relpaths[prediction_id]
        break ### TMP

    settings['rosetta_args_list_list'][len(job_dicts)].extend(
        repack_generic_args
    )
    settings['rosetta_args_list_list'][len(job_dicts)].extend( [
        '-out:prefix mutate_',
    ] )
    settings['appname_list'].append( 'fixbb' )
    job_dicts.append( job_dict )

    # Minimize repacked wildtype step
    job_dict = {}
    for prediction_id in original_job_dict:
        data_dir = os.path.join( os.path.join(output_dir, 'data'), '%d' % prediction_id )
        pdb_id = pdb_ids_from_pred_ids[prediction_id]
        for n_struct in xrange(1, inner_jobs + 1):
            inner_prediction_id = '%s/%d-%d-minimize-repack_wt' % (str(prediction_id), n_struct, len(job_dicts))
            repacked_inner_prediction_id = '%s/%d-%d-repack_wt' % (str(prediction_id), n_struct, len(job_dicts)-2)

            inner_prediction_repacked_dir = os.path.join(output_dir, repacked_inner_prediction_id)
            repacked_pdb = os.path.relpath(os.path.join(inner_prediction_repacked_dir, 'repack-wt_bkrb_min_cst.%s_0001_0001_last_0001.pdb.gz' % pdb_id), output_dir)
            job_dict[inner_prediction_id] = {}
            job_dict[inner_prediction_id]['input_file_list'] = [
                repacked_pdb
            ]
        break ### TMP

    settings['rosetta_args_list_list'][len(job_dicts)].extend(
        generic_min_args
    )
    settings['appname_list'].append( 'minimize_with_cst' )
    job_dicts.append( job_dict )

    # Minimize repacked mutant step
    job_dict = {}
    for prediction_id in original_job_dict:
        data_dir = os.path.join( os.path.join(output_dir, 'data'), '%d' % prediction_id )
        pdb_id = pdb_ids_from_pred_ids[prediction_id]
        for n_struct in xrange(1, inner_jobs + 1):
            inner_prediction_id = '%s/%d-%d-minimize-mutate' % (str(prediction_id), n_struct, len(job_dicts))
            repacked_inner_prediction_id = '%s/%d-%d-mutate' % (str(prediction_id), n_struct, len(job_dicts)-2)

            inner_prediction_repacked_dir = os.path.join(output_dir, repacked_inner_prediction_id)
            repacked_pdb = os.path.relpath(os.path.join(inner_prediction_repacked_dir, 'mutate_bkrb_min_cst.%s_0001_0001_last_0001.pdb.gz' % pdb_id), output_dir)
            job_dict[inner_prediction_id] = {}
            job_dict[inner_prediction_id]['input_file_list'] = [
                repacked_pdb
            ]
        break ### TMP

    settings['rosetta_args_list_list'][len(job_dicts)].extend(
        generic_min_args
    )
    settings['appname_list'].append( 'minimize_with_cst' )
    job_dicts.append( job_dict )

    # Rescore minimized wt by separating partners
    # job_dict = {}
    # chains_to_move_cache = {}
    # score_fxn = 'talaris2013' ### Note WARNING
    # for prediction_id in original_job_dict:
    #     pdb_id = pdb_ids_from_pred_ids[prediction_id]
    #     job_details = ppi_api.get_job_details(prediction_id)
    #     substitution_parameters = json.loads(job_details['JSONParameters'])
    #     chains_to_move = substitution_parameters['%%chainstomove%%']
    #     chains_to_move_cache[prediction_id] = chains_to_move
    #     for n_struct in xrange(1, inner_jobs + 1):
    #         inner_prediction_id = '%s/%d-%d-rescore_sep_wt' % (str(prediction_id), n_struct, len(job_dicts))
    #         finalminned_inner_prediction_id = '%s/%d-%d-minimize-repack_wt' % (str(prediction_id), n_struct, len(job_dicts)-2)

    #         inner_prediction_finalminned_dir = os.path.join(output_dir, finalminned_inner_prediction_id)
    #         finalminned_pdb = os.path.relpath(os.path.join(inner_prediction_finalminned_dir, 'mutate_bkrb_min_cst.%s_0001_0001_last_0001.pdb.gz' % pdb_id), output_dir)
    #         job_dict[inner_prediction_id] = {}
    #         job_dict[inner_prediction_id]['input_file_list'] = [
    #             finalminned_pdb
    #         ]
    #     break ### TMP

    # settings['rosetta_args_list_list'][len(job_dicts)].extend(
    #     generic_min_args
    # )
    # settings['appname_list'].append( 'minimize_with_cst' )
    # job_dicts.append( job_dict )

    print job_dicts ### TMP

    assert( len(job_dicts) == num_steps )
    ct = cluster_template.ClusterTemplate(num_steps, settings_dict = settings)
    for i, job_dict in enumerate(job_dicts):
        ct.set_job_dict(job_dict, step_num = i)

    ct.write_runs()
    print 'Job files written to directory:', os.path.abspath(output_dir)
