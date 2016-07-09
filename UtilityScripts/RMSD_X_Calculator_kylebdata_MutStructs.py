# This script is intended to calculate the average RMSD of each Rosetta-generated mutant ensemble member to the known
# mutant crystal structure

from Bio.PDB import *
import pyRMSD
import pyRMSD.RMSDCalculator
from pyRMSD.availableCalculators import availableCalculators
import scipy.spatial.distance
import prody
import numpy as np
import pprint
import sys
import re
import pandas as pd
import os
import tempfile
import zipfile
import shutil
import multiprocessing
from klab.bio.clustalo import PDBSeqresSequenceAligner

from kddg.api import settings
sys_settings = settings.load()

# From Kyle's Finalize.py
def find_neighbors(mutations, open_strct, neighbor_distance=8.0):
    # mutations = mutations_resfile(filenum_dir)
    # There should only be one model in PDB file
    num_models = 0
    for model in open_strct.get_models():
        num_models += 1
    assert (num_models == 1)

    chain_list = [chain.get_id() for chain in open_strct[0].get_chains()]
    neighbors_set = set()
    for mutation in mutations:
        res_id, chain_id, pikaa, mut_aa = mutation
        mut_chain = str(chain_id)
        try:
            mut_pos = int(res_id)
            mut_insertion_code = ' '
        except ValueError:
            mut_pos = int(res_id[:-1])
            mut_insertion_code = res_id[-1]

        mut_residue = open_strct[0][mut_chain][(' ', mut_pos, mut_insertion_code)]
        for chain in chain_list:
            for residue in [res.get_id() for res in open_strct[0][chain].get_residues()]:
                try:
                    # Kyle note - might be good to do something else for consistency, since not all residues have CB
                    dist = mut_residue['CB'] - open_strct[0][chain][residue]['CB']
                    if dist < neighbor_distance:
                        neighbors_set.add((residue, chain))
                except KeyError:
                    try:
                        dist = mut_residue['CA'] - open_strct[0][chain][residue]['CA']
                        if dist < neighbor_distance:
                            neighbors_set.add((residue, chain))
                    except KeyError:
                        pass
    neighbors = [set_residue for set_residue in neighbors_set]

    return neighbors

# Kyle's RMSDWrapper.py (Spread out all over the place now)
def rmsd(pyrmsd_calc, coordinates, fresh_coords):
    def generate_out_dict(rmsd):
        # from pyRMSD.condensedMatrix import CondensedMatrix
        # out_dict_temp['Mean'] = CondensedMatrix(rmsd).calculateMean()
        # out_dict_temp['Variance'] = CondensedMatrix(rmsd).calculateVariance()
        # out_dict_temp['Skewness'] = CondensedMatrix(rmsd).calculateSkewness()
        # out_dict_temp['Kurtosis'] = CondensedMatrix(rmsd).calculateKurtosis()
        # out_dict_temp['Raw'] = CondensedMatrix(rmsd).get_data().flatten()

        out_dict_temp = {}
        out_dict_temp['Mean'] = np.mean(rmsd)
        out_dict_temp['Variance'] = np.var(rmsd)
        out_dict_temp['Raw'] = np.asarray(rmsd)
        return out_dict_temp

    def superpose(coordinates, fresh_coords):
        # Superimposes input_PDBs onto reference structures using only backbone atoms
        # Select reference PDB backbone atoms
        ref_coords = fresh_coords[0].select('name CA C N').getCoords()

        # http://nghiaho.com/?page_id=671
        # http://nghiaho.com/uploads/code/rigid_transform_3D.py_
        coordinates_superposed_tmp = []
        for coordinate in coordinates:
            fitting_coords = np.asarray(coordinate.select('name CA C N').getCoords())

            centroid_fitting = np.mean(fitting_coords, axis=0)
            centroid_ref = np.mean(ref_coords, axis=0)

            centered_fitting = fitting_coords - np.tile(centroid_fitting, (fitting_coords.shape[0], 1))
            centered_ref = ref_coords - np.tile(centroid_ref, (ref_coords.shape[0], 1))

            H = np.matmul(centered_fitting.T, centered_ref)
            U, S, Vt = np.linalg.svd(H)

            # List containing Rotation Matrix
            Rotation_matrix = np.matmul(Vt.T, U.T)
            Translation_matrix = -np.matmul(Rotation_matrix, centroid_fitting.T) + centroid_ref.T

            # CHECK THIS!!! Deviated from reference script
            coordinate_superposed = (np.matmul(Rotation_matrix, coordinate.getCoords().T) + np.tile(Translation_matrix,(coordinate.getCoords().shape[0], 1)).T).T
            coordinates_superposed_tmp.append(coordinate_superposed)

        # # Double-checking alignment with a 3D plot
        # import matplotlib.pyplot as plt
        # from mpl_toolkits.mplot3d import Axes3D
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # for coordinate_superposed in coordinates_superposed_tmp:
        #     x = coordinate_superposed[:, 0]
        #     y = coordinate_superposed[:, 1]
        #     z = coordinate_superposed[:, 2]
        #     ax.scatter(x, y, -z, zdir='z', c='red')
        # a = ref_coords[:, 0]
        # b = ref_coords[:, 1]
        # c = ref_coords[:, 2]
        # ax.scatter(a, b, -c, zdir='z', c='blue', s=100)
        # plt.show()

        return np.asarray(coordinates_superposed_tmp)

    if type(coordinates) is dict:
        rmsd_list = {}
        for mutres, coord_set in coordinates.iteritems():
            coordinates_superposed = superpose(coord_set, fresh_coords[mutres])
            rmsd_array = []
            for rosettaout in coordinates_superposed:
                error = rosettaout - np.asarray([fresh_coords[mutres][0].getCoords()])
                rmsd_array.append(np.sqrt(np.sum(np.multiply(error, error)) / fresh_coords[mutres][0].getCoords().shape[0]))

            out_dict = generate_out_dict(rmsd_array)
            rmsd_list[mutres] = out_dict

            # coords_plus_ref = np.concatenate((coord_set, fresh_coords[mutres]), axis=0)
            # calculator = pyRMSD.RMSDCalculator.RMSDCalculator(pyrmsd_calc, coords_plus_ref)
            # rmsd = calculator.oneVsTheOthers(50)
            # out_dict = generate_out_dict(rmsd)
            # rmsd_list.append([mutres, out_dict])

        return rmsd_list
    else:

        # # http://prody.csb.pitt.edu/tutorials/ensemble_analysis/dimer.html#multimeric-structures
        # # Set up Ensemble
        # ensemble = prody.PDBEnsemble('ensemble')
        # ensemble.setAtoms(fresh_coords_backbone)
        # ensemble.setCoords(fresh_coords_backbone.getCoords())
        # for coordinate in coordinates:
        #     ensemble.addCoordset(coordinate)
        # ensemble.iterpose()
        # # Convert superposed PDBs to numpy array coordinates
        # superposed_PDBs = list([])
        # for superposed_PDB in ensemble:
        #     superposed_PDBs.append(superposed_PDB.getCoords())
        # coordinates_np = np.asarray(superposed_PDBs)

        coordinates_superposed = superpose(coordinates, fresh_coords)

        rmsd_array = []
        for rosettaout in coordinates_superposed:
            error = rosettaout - np.asarray([fresh_coords[0].getCoords()])
            rmsd_array.append(np.sqrt(np.sum(np.multiply(error, error))/fresh_coords[0].getCoords().shape[0]))

        # Old stuff
        # coordinates_np = np.asarray([coordinate.getCoords() for coordinate in coordinates])
        # coords_plus_ref = np.concatenate((coordinates_np, np.asarray([fresh_coords[0].getCoords()])), axis=0)

        # pyRMSD Method of calculating RMSDs - [50] is relevant row of the RMSD Matrix
        # coords_plus_ref = np.concatenate((coordinates_superposed, np.asarray([fresh_coords[0].getCoords()])), axis=0)
        # calculator = pyRMSD.RMSDCalculator.RMSDCalculator(pyrmsd_calc, coords_plus_ref)
        # rmsd = calculator.pairwiseRMSDMatrix()
        # out_dict = generate_out_dict(rmsd)

        out_dict = generate_out_dict(rmsd_array)
        return out_dict

def flip_me(residue_maps, wt_to_mut_chains):
    # Flips residue_maps and wt_to_mut_chains for wt --> mut residue numbering
    residue_maps_reverse = {}
    for nested_dict in residue_maps:
        residue_maps_reverse[(nested_dict[1], nested_dict[0])] = {mut: wt for wt, mut in
                                                                  residue_maps[nested_dict].iteritems()}
    # Flips residue_maps and wt_to_mut_chains for wt --> mut residue numbering
    mut_to_wt_chains = {mut: wt for wt, mut in wt_to_mut_chains.iteritems()}
    return residue_maps_reverse, mut_to_wt_chains

def common_atoms(tmp_mut_pdb, tmp_wt_pdb, residue_maps, wt_to_mut_chains, input_pdbs, selection_type):
    # Checks for atoms present in both RosettaOut and Reference mutant structures
    mut_atoms = set()
    wt_atoms = set()
    real_wt_atoms = set()
    real_mut_atoms = set()
    acceptable_atoms_wt_set = set()
    acceptable_atoms_mut_set = set()

    if selection_type == 'CA':
        for chain_mapping in residue_maps:
            for wt_residue, mut_residue in residue_maps[chain_mapping].iteritems():
                wt_atoms.add(wt_residue.split()[0] + wt_residue.split()[1])
                mut_atoms.add(mut_residue.split()[0] + mut_residue.split()[1])

        mut_pdb_hv = prody.parsePDB(tmp_mut_pdb).select('calpha').getHierView()
        wt_pdb_hv = prody.parsePDB(input_pdbs[1]).select('calpha').getHierView()

        # Intersection between Mut dict atoms and Mut hv atoms = mut real atoms
        for hv_chain in mut_pdb_hv:
            for hv_calpha in hv_chain:
                for (chain, position) in zip(hv_calpha.getChids(), hv_calpha.getResnums()):
                    if (chain + str(position)) in mut_atoms:
                        real_mut_atoms.add(chain + str(position))
        # Intersection between WT dict atoms and WT hv atoms = wt real atoms
        for hv_chain in wt_pdb_hv:
            for hv_calpha in hv_chain:
                for (chain, position) in zip(hv_calpha.getChids(), hv_calpha.getResnums()):
                    if (chain + str(position)) in wt_atoms:
                        real_wt_atoms.add(chain + str(position))

        # Intersection between wt and mut real atoms (use wt_to_mut_dict)
        for real_wt_atom in real_wt_atoms:
            chain = real_wt_atom[:1]
            position = real_wt_atom[1:]
            if (wt_to_mut_chains[chain] + str(residue_maps[(chain, wt_to_mut_chains[chain])]['%s %s ' % (chain, ('   ' + str(position))[-3:])].split()[1])) in real_mut_atoms:
                acceptable_atoms_wt_set.add(real_wt_atom)
                acceptable_atoms_mut_set.add(wt_to_mut_chains[chain] + str(residue_maps[(chain, wt_to_mut_chains[chain])]['%s %s ' % (chain, ('   ' + str(position))[-3:])].split()[1]))

        acceptable_atoms = {}
        acceptable_atoms['RosettaOut'] = acceptable_atoms_wt_set
        acceptable_atoms['Mutant PDB'] = acceptable_atoms_mut_set

        return acceptable_atoms

    # For point mutants and neighbors (only looks at residues returned by find_neighbors(), doesn't look at +1/-1 residues)
    if isinstance(selection_type, list):
        mut_pdb_hv = prody.parsePDB(tmp_mut_pdb).getHierView()
        wt_pdb_hv = prody.parsePDB(input_pdbs[1]).getHierView()

        for residue in selection_type:
            wt_res = wt_pdb_hv[(residue[1], int(residue[0]))]
            try:
                mut_res = mut_pdb_hv[(wt_to_mut_chains[residue[1]], int(residue_maps[(residue[1], wt_to_mut_chains[residue[1]])]['%s %s ' % (residue[1], ('   ' + str(residue[0]))[-3:])].split()[1]))]
                if wt_res.getResname() == mut_res.getResname():
                    for wt_atom in wt_res:
                        for mut_atom in mut_res:
                            if wt_atom.getName() == mut_atom.getName():
                                acceptable_atoms_wt_set.add((wt_res.getChid(), wt_res.getResname(), wt_res.getResnum(), wt_atom.getName()))
                                acceptable_atoms_mut_set.add((mut_res.getChid(), mut_res.getResname(), mut_res.getResnum(), mut_atom.getName()))
            except:
                print '%s%s has no equivalent residue in the Reference Mutant PDB' %(residue[1], int(residue[0]))
                pass

        return acceptable_atoms_wt_set, acceptable_atoms_mut_set

def global_ca_coordinates(input_pdbs, tmp_mut_pdb, tmp_wt_pdb, residue_maps, wt_to_mut_chains, input_type):
    acceptable_atoms = common_atoms(tmp_mut_pdb, tmp_wt_pdb, residue_maps, wt_to_mut_chains, input_pdbs, selection_type='CA')

    def generate_ca_atom_list(input_pdbs, acceptable_atoms, input_type):
        temp = []
        # Adds acceptable atom indexes to atom list which are then used in select() to return atom object with desired coordinates
        for input_pdb in input_pdbs:
            if 'WT.' not in input_pdb:
                atom_list = []
                c_alpha_all = prody.parsePDB(input_pdb).select('calpha')
                c_alpha_all_hv = c_alpha_all.getHierView()
                for c_alpha_chain in c_alpha_all_hv:
                    for c_alpha in c_alpha_chain:
                        if input_type == 'Mutant PDB':
                            if c_alpha.getChids()[0] + str(c_alpha.getResnums()[0]) in acceptable_atoms['Mutant PDB']:
                                for atom in c_alpha:
                                    atom_list.append(str(atom.getIndex()))
                        if input_type == 'RosettaOut':
                            if c_alpha.getChids()[0] + str(c_alpha.getResnums()[0]) in acceptable_atoms['RosettaOut']:
                                for atom in c_alpha:
                                    atom_list.append(str(atom.getIndex()))
                temp.append(c_alpha_all.select('index ' + ' '.join(atom_list)))

        return temp

    if input_type == 'RosettaOut':
        return generate_ca_atom_list(input_pdbs, acceptable_atoms, input_type), acceptable_atoms
    if input_type == 'Mutant PDB':
        return generate_ca_atom_list([tmp_mut_pdb], acceptable_atoms, input_type)

def neighborhood_coordinates(neighbors, input_pdbs, residue_maps, wt_to_mut_chains, tmp_mut_pdb, tmp_wt_pdb, input_type):
    neighbors_reorganized = [[neighbor[0][1], neighbor[1]] for neighbor in neighbors]
    acceptable_atoms_wt_set, acceptable_atoms_mut_set = common_atoms(tmp_mut_pdb, tmp_wt_pdb, residue_maps, wt_to_mut_chains, input_pdbs, neighbors_reorganized)
    residue_maps_reverse, mut_to_wt_chains = flip_me(residue_maps, wt_to_mut_chains)

    # Converts neighbor selection into mutant PDB numbering
    if input_type == 'Mutant PDB':
        neighbors_mut_numbering = []

        for neighbor_residue in neighbors:
            try:
                neighbors_mut_numbering.append(((' ',
                                         int(residue_maps[(neighbor_residue[1], wt_to_mut_chains[neighbor_residue[1]])][
                                                 '%s %s ' % (
                                                 neighbor_residue[1], ('   ' + str(neighbor_residue[0][1]))[-3:])].split()[1]),
                                         ' '),
                                        residue_maps[(neighbor_residue[1], wt_to_mut_chains[neighbor_residue[1]])]['%s %s ' % (
                                        neighbor_residue[1], ('   ' + str(neighbor_residue[0][1]))[-3:])].split()[0]))
            except:
                print '%s%s has no equivalent residue in the Reference Mutant PDB' %(neighbor_residue[1], str(neighbor_residue[0][1]))
                pass

    def generate_neighborhood_atom_list(input_pdbs, neighbors, acceptable_atoms_wt_set, acceptable_atoms_mut_set, input_type):
        coordinates = []
        for input_pdb in input_pdbs:
            if 'WT.' not in input_pdb:
                atom_list = []
                neighborhood = prody.parsePDB(input_pdb)
                neighborhood_hv = neighborhood.getHierView()
                res_list = [neighborhood_hv[neighbor[1], neighbor[0][1]] for neighbor in neighbors]

                for res in res_list:
                    # Check if numbering should be for WT or Mutant
                    # Check that residues are present in acceptable_residues
                    # Check that atom coordinates are present in acceptable_atoms
                    for atom in res:
                        if input_type == 'Mutant PDB':
                            # print (mut_to_wt_chains[res.getChid()], res.getResname(),
                            #     int(residue_maps_reverse[(res.getChid(), mut_to_wt_chains[res.getChid()])][
                            #             '%s %s ' % (res.getChid(), ('   ' + str(res.getResnum()))[-3:])].split()[1]),
                            #     atom.getName())
                            if (res.getChid(), res.getResname(),int(res.getResnum()), atom.getName()) in acceptable_atoms_mut_set:
                                if atom.getElement() != 'H':
                                    atom_list.append(str(atom.getIndex()))

                        else:
                            if (res.getChid(), res.getResname(), res.getResnum(), atom.getName()) in acceptable_atoms_wt_set:
                                if atom.getElement() != 'H':
                                    atom_list.append(str(atom.getIndex()))

                coordinates.append(neighborhood.select('index ' + ' '.join(atom_list)))

        return coordinates

    if input_type == 'RosettaOut':
        return generate_neighborhood_atom_list(input_pdbs, neighbors,  acceptable_atoms_wt_set, acceptable_atoms_mut_set, input_type), acceptable_atoms_wt_set
    if input_type == 'Mutant PDB':
        return generate_neighborhood_atom_list([tmp_mut_pdb], neighbors_mut_numbering,  acceptable_atoms_wt_set, acceptable_atoms_mut_set, input_type)

def mutant_coordinates(input_pdbs, mutations, residue_maps, wt_to_mut_chains, tmp_mut_pdb, tmp_wt_pdb, input_type):
    acceptable_atoms_wt_set, acceptable_atoms_mut_set = common_atoms(tmp_mut_pdb, tmp_wt_pdb, residue_maps, wt_to_mut_chains, input_pdbs, mutations)
    mut_key_dict = {}
    if input_type == 'Mutant PDB':
        mutations_mut_numbering = []
        for mutation in mutations:
        # Converts mutation in WT PDB numbering to Mutant PDB numbering
            fetch_from_res_map = residue_maps[(mutation[1], wt_to_mut_chains[mutation[1]])]['%s %s ' % (mutation[1], ('   ' + str(mutation[0]))[-3:])].split()
            mutations_mut_numbering.append([fetch_from_res_map[1], fetch_from_res_map[0]])
            mut_key_dict[str(fetch_from_res_map[0]) + fetch_from_res_map[1]] = mutation[1] + str(mutation[0])

    def generate_point_atom_list(input_pdbs, mutations, acceptable_atoms_wt_set, acceptable_atoms_mut_set, mut_key_dict, input_type):
        mutation_dict = {}
        for counter, mutation in enumerate(mutations):
            temp = []
            for input_pdb in input_pdbs:
                if 'WT.' not in input_pdb:
                    atom_list = []
                    point_mutant = prody.parsePDB(input_pdb)
                    point_mutant_hv = point_mutant.getHierView()
                    res_list = [point_mutant_hv[mutation[1], int(mutation[0])]]
                    for res in res_list:
                        # Check if numbering should be for WT or Mutant
                        # Check that residues are present in acceptable_residues
                        # Check that atom coordinates are present in acceptable_atoms
                        for atom in res:
                            if input_type == 'Mutant PDB':
                                if (res.getChid(), res.getResname(),int(res.getResnum()), atom.getName()) in acceptable_atoms_mut_set:
                                    if atom.getElement() != 'H':
                                        atom_list.append(str(atom.getIndex()))
                            if input_type == 'RosettaOut':
                                if (res.getChid(), res.getResname(), res.getResnum(), atom.getName()) in acceptable_atoms_wt_set:
                                    if atom.getElement() != 'H':
                                        atom_list.append(str(atom.getIndex()))

                    temp.append(point_mutant.select('index ' + ' '.join(atom_list)))
            if input_type == 'Mutant PDB':
                mutation_dict[mut_key_dict[mutation[1] + str(mutation[0])]] = temp
            if input_type == 'RosettaOut':
                mutation_dict[mutation[1] + str(mutation[0])] = temp
        return mutation_dict

    if input_type == 'RosettaOut':
        return generate_point_atom_list(input_pdbs, mutations, acceptable_atoms_wt_set, acceptable_atoms_mut_set, mut_key_dict, input_type), acceptable_atoms_wt_set
    if input_type == 'Mutant PDB':
        return generate_point_atom_list([tmp_mut_pdb], mutations_mut_numbering, acceptable_atoms_wt_set, acceptable_atoms_mut_set, mut_key_dict, input_type)

def bin_me(templist):
    from collections import Counter
    bins = {}
    floor_list = [X // 30 for X in templist]
    X_counts = Counter(floor_list)
    bin_ranges = {int(-6.0): '-180 <= x < -150',
                  int(-5.0): '-150 <= x < -120',
                  int(-4.0): '-120 <= x < -90',
                  int(-3.0): '-90 <= x < -60',
                  int(-2.0): '-60 <= x <-30',
                  int(-1.0): '-30 <= x < 0',
                  int(0): '0 <= x < 30',
                  int(1.0): '30 <= x < 60',
                  int(2.0): '60 <= x < 90',
                  int(3.0): '90 <= x < 120',
                  int(4.0): '120 <= x < 150',
                  int(5.0): '150 <= x <180',
                  int(6.0): 'x = 180'}
    readable_X_counts = dict((bin_ranges[key], value) for (key, value) in X_counts.items())
    return readable_X_counts

def PCA_Analysis(predID, input_pdbs, tmp_mut_pdb):
    # http://prody.csb.pitt.edu/tutorials/ensemble_analysis/xray_calculations.html#calculations
    ref_structure = prody.parsePDB(tmp_mut_pdb)
    ref_selection = ref_structure.select('calpha')
    for ref_chain in ref_structure.getHierView():
        prody.startLogfile('%s_%s_log'%(predID, ref_chain.getChid()))
        ensemble = prody.PDBEnsemble('%s_%s_ensemble'%(predID, ref_chain.getChid()))
        ensemble.setAtoms(ref_chain.select('calpha'))
        ensemble.setCoords(ref_chain.select('calpha'))
        print 'Generating Ensemble...'
        for input_pdb in input_pdbs:
            structure = prody.parsePDB(input_pdb, subset='calpha')
            mappings = prody.mapOntoChain(structure, ref_chain)
            atommap = mappings[0][0]
            ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
        ensemble.iterpose()
        prody.closeLogfile('%s_%s_log'%(predID, ref_chain.getChid()))
        prody.writePDB('%s_%s_ensemble.pdb'%(predID, ref_chain.getChid()), ensemble)

        pca = prody.PCA('%s_%s_pca'%(predID, ref_chain.getChid()))
        pca.buildCovariance(ensemble)
        pca.calcModes()
        print pca.getEigvecs()

def chi_angles(input_pdbs, mutations):
    # http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
    chi1_dict = {'ARG': ['N', 'CA', 'CB', 'CG'],
                 'ASN': ['N', 'CA', 'CB', 'CG'],
                 'ASP': ['N', 'CA', 'CB', 'CG'],
                 'CYS': ['N', 'CA', 'CB', 'SG'],
                 'GLN': ['N', 'CA', 'CB', 'CG'],
                 'GLU': ['N', 'CA', 'CB', 'CG'],
                 'HIS': ['N', 'CA', 'CB', 'CG'],
                 'ILE': ['N', 'CA', 'CB', 'CG1'],
                 'LEU': ['N', 'CA', 'CB', 'CG'],
                 'LYS': ['N', 'CA', 'CB', 'CG'],
                 'MET': ['N', 'CA', 'CB', 'CG'],
                 'PHE': ['N', 'CA', 'CB', 'CG'],
                 'PRO': ['N', 'CA', 'CB', 'CG'],
                 'SER': ['N', 'CA', 'CB', 'OG'],
                 'THR': ['N', 'CA', 'CB', 'OG1'],
                 'TRP': ['N', 'CA', 'CB', 'CG'],
                 'TYR': ['N', 'CA', 'CB', 'CG'],
                 'VAL': ['N', 'CA', 'CB', 'CG1']}
    chi2_dict = {'ARG': ['CA', 'CB', 'CG', 'CD'],
                 'ASN': ['CA', 'CB', 'CG', 'OD1'],
                 'ASP': ['CA', 'CB', 'CG', 'OD1'],
                 'GLN': ['CA', 'CB', 'CG', 'CD'],
                 'GLU': ['CA', 'CB', 'CG', 'CD'],
                 'HIS': ['CA', 'CB', 'CG', 'ND1'],
                 'ILE': ['CA', 'CB', 'CG1', 'CD1'],
                 'LEU': ['CA', 'CB', 'CG', 'CD1'],
                 'LYS': ['CA', 'CB', 'CG', 'CD'],
                 'MET': ['CA', 'CB', 'CG', 'SD'],
                 'PHE': ['CA', 'CB', 'CG', 'CD1'],
                 'PRO': ['CA', 'CB', 'CG', 'CD'],
                 'TRP': ['CA', 'CB', 'CG', 'CD1'],
                 'TYR': ['CA', 'CB', 'CG', 'CD1']}

    xangles_dict = {}

    for mutation in mutations:
        if mutation[3] == 'A' or mutation[3] == 'G':
            xangles_dict['%s%s' % (mutation[1], mutation[0])] = "Mutation is %s!" % mutation[3]
        else:
            xangles_dict['%s%s' % (mutation[1], mutation[0])] = {}
            x1_templist = []
            x2_templist = []
            for input_pdb in input_pdbs:
                if 'WT.pdb_' not in input_pdb:
                    p = prody.parsePDB(input_pdb)
                    hv = p.getHierView()

                    current_res = hv[mutation[1], int(mutation[0])]

                    if current_res.getResname() in chi1_dict.keys():
                        current_res.select('name %s' % chi1_dict[current_res.getResname()])
                        chi1 = prody.measure.measure.calcDihedral(
                            current_res.select('name %s' % chi1_dict[current_res.getResname()][0]),
                            current_res.select('name %s' % chi1_dict[current_res.getResname()][1]),
                            current_res.select('name %s' % chi1_dict[current_res.getResname()][2]),
                            current_res.select('name %s' % chi1_dict[current_res.getResname()][3]))
                        x1_templist.append(chi1[0])

                    if current_res.getResname() in chi2_dict.keys():
                        chi2 = prody.measure.measure.calcDihedral(
                            current_res.select('name %s' % chi2_dict[current_res.getResname()][0]),
                            current_res.select('name %s' % chi2_dict[current_res.getResname()][1]),
                            current_res.select('name %s' % chi2_dict[current_res.getResname()][2]),
                            current_res.select('name %s' % chi2_dict[current_res.getResname()][3]))
                        x2_templist.append(chi2[0])

            # Adds raw X angles to output dictionary
            xangles_dict['%s%s' % (mutation[1], mutation[0])]['X1'] = x1_templist

            if x2_templist != []:
                xangles_dict['%s%s' % (mutation[1], mutation[0])]['X2'] = x2_templist

            # x1_bins = bin_me(x1_templist)
            # xangles_dict['%s%s' % (mutation[0], mutation[1])]['X1'] = x1_bins
            #
            # if x2_templist != []:
            #     x2_bins = bin_me(x2_templist)
            #     xangles_dict['%s%s' % (mutation[0], mutation[1])]['X2'] = x2_bins

                # Outputs Ramachandran-like plots for X1 and X2, Violin plots (in progress) if only X1 is present
                # import pandas as pd
                # import seaborn as sns
                # import matplotlib.pyplot as plt

                # x1 = pd.Series(x1_templist, name = 'X1')
                # if x2_templist != []:
                #    x2 = pd.Series(x2_templist, name = 'X2')
                #    mrplotty = sns.regplot(x1, x2, fit_reg = False)
                #    plt.savefig('%s%splot.pdf'  %(mutation[0], mutation[1]))

    return xangles_dict

# Return mutations for Prediction ID

def Fetch_Mutant_ID(predID):
    # List of dictionaries for Resfile mutation info
    sys.path.insert(0,'/kortemmelab/home/james.lucas')  # this should point to the directory above your ddg repo checkout
    from kddg.api.ppi import get_interface_with_config_file as get_ppi_interface

    # Create an interface to the database
    ppi_api = get_ppi_interface()

    # Get details back for one prediction
    PredID_Details = ppi_api.get_job_details(predID, include_files=True)

    mutations = []
    for mutation_entry in PredID_Details['PDBMutations']:
        mutations.append([mutation_entry['ResidueID'].strip(),
                          mutation_entry['Chain'],
                          mutation_entry['PDBFileID'],
                          mutation_entry['MutantAA']])

    # mutant_list contains (Wildtype_PDB_ID, MUTANT_PDB_ID, MutagenesisID_index)
    MutagenesisID_list = []
    df = pd.read_csv('/kortemmelab/home/james.lucas/skempi_mutants.tsv', delimiter='\t')
    for index, row in df.iterrows():
        if row['PPMutagenesisID'] == PredID_Details['PDBMutations'][0]['PPMutagenesisID']:
            MutagenesisID_list.append((row['Wildtype'], row['Mutant'], index))
            break

    #Grabs Rosetta energies for complexes
    scores = ppi_api.get_prediction_scores(predID)
    for score_method in scores:
        REU_list = [scores[score_method][ensemble_member]['MutantComplex']['total'] for ensemble_member in scores[score_method]]

    return MutagenesisID_list, mutations, PredID_Details, df, REU_list

def Generate_PDBs_and_Resmaps(PredID_Details, Mutant_PDB_ID, MutagenesisID_index, df):
    wt_pdb_filename = PredID_Details['Files']['Input'][0]['Filename']  # Literally prints the filename (1TM1_EI.pdb)
    wt_pdb_chains = [chain for chain in re.sub('_|\.', ' ', wt_pdb_filename).split()[1]]  # Prints chains as string (EI)

    # WT pdb is from a string and Mut PDB is from file... so... yeah... this happened
    raw_wt_pdb = PredID_Details['Files']['Input'][0]['Content']
    raw_mutant_pdb = ''
    for line in open('/kortemmelab/home/james.lucas/skempi_mutants/%s.pdb' % Mutant_PDB_ID):
        raw_mutant_pdb += line

    WT_PDB_ID = wt_pdb_filename[:4]

    # IN PROGRESS
    residue_maps, wt_to_mut_chains = map_pdbs(WT_PDB_ID, Mutant_PDB_ID, wt_pdb_chains, df, MutagenesisID_index)  # Variable df is temporary while Shane figures out chain mapping thing
    fresh_wt_pdb, tmp_wt_pdb = strip_pdbs(raw_wt_pdb, wt_to_mut_chains, input_type='WT PDB')
    fresh_mut_pdb, tmp_mut_pdb = strip_pdbs(raw_mutant_pdb, wt_to_mut_chains, input_type='Mutant PDB')

    return fresh_wt_pdb, tmp_wt_pdb, fresh_mut_pdb, tmp_mut_pdb, residue_maps, wt_to_mut_chains

def strip_pdbs(raw_pdb, wt_to_mut_chains, input_type):
    import tempfile

    if input_type == 'Mutant PDB':
        pdb_chains = []
        for key, item in wt_to_mut_chains.items():
            pdb_chains.append(item)

    else:
        pdb_chains = []
        for key, item in wt_to_mut_chains.items():
            pdb_chains.append(key)

    with tempfile.NamedTemporaryFile(mode='w+b', delete=False) as temp_PDBFile:
        tmp_pdb = temp_PDBFile.name
        # for line in raw_pdb.splitlines(): #Switched from string to file, no longer need to split string for lines
        # lol nvm
        for line in raw_pdb.splitlines():
            parsed = line.split()
            if parsed[0].strip() == 'ATOM':
                if parsed[4] in pdb_chains:
                    temp_PDBFile.write(line + '\n')
            elif parsed[0].strip() == 'HETATM':
                if parsed[4] in pdb_chains:
                    temp_PDBFile.write(line + '\n')
            else:
                temp_PDBFile.write(line + '\n')

    with open(tmp_pdb, mode='rb+') as temp_PDBFile:
        # WARNING: BIOPYTHON!!!
        parser = PDBParser(PERMISSIVE=1)
        fresh_pdb = parser.get_structure('Open', temp_PDBFile.name)

    return fresh_pdb, tmp_pdb

def map_pdbs(WT_PDB_ID, Mutant_PDB_ID, wt_pdb_chains, df, MutagenesisID_index):
    pssa = PDBSeqresSequenceAligner(WT_PDB_ID, Mutant_PDB_ID, cache_dir='/tmp')
    mut_chain_ids = {}
    residue_maps = {}

    for wt_chain_id in wt_pdb_chains:
        mut_chain_ids[wt_chain_id] = pssa.get_matching_chains(wt_chain_id)
        for mut_chain_id in mut_chain_ids[wt_chain_id]:
            residue_maps[(wt_chain_id, mut_chain_id)] = pssa.get_atom_residue_mapping(wt_chain_id, mut_chain_id)

    # Tempory chain matching solution...
    wt_to_mut_chains = {}
    for mapping in re.sub(' |,', ' ', df.loc[MutagenesisID_index]['WT:Mut Mapping']).split():
        wt_to_mut_chains[mapping[0]] = mapping[2]

    chains_to_keep = set([key for key in residue_maps.keys()]) & set([(wt_to_mut , wt_to_mut_chains[wt_to_mut]) for wt_to_mut in wt_to_mut_chains])

    for residue_maps_key in residue_maps.keys():
        if residue_maps_key not in chains_to_keep:
            residue_maps.pop(residue_maps_key)

    return residue_maps, wt_to_mut_chains

def do_math(outputdir, predID):
    # # Use CUDA for GPU calculations, if avialable
    # if 'QCP_CUDA_MEM_CALCULATOR' in availableCalculators():
    #     pyrmsd_calc = 'QCP_CUDA_MEM_CALCULATOR'
    # else:
    #     pyrmsd_calc = 'QCP_SERIAL_CALCULATOR'

    # pyRMSD calculator which does not perform superposition
    pyrmsd_calc = 'NOSUP_SERIAL_CALCULATOR'

    input_temp = []

    for app_output_dir in os.listdir(outputdir):
        for app_outfile in os.listdir(os.path.join(outputdir, app_output_dir)):
            # Mutant Structures
            if app_outfile == 'min_cst.mutate_bkrb_min_cst_0001_0001.pdb.gz':
                input_temp.append(os.path.join(outputdir, app_output_dir, app_outfile))
                input_pdbs = sorted(input_temp)
                # WT Structures
                # if app_outfile == 'min_cst.repack-wt_bkrb_min_cst_0001_0001.pdb.gz':
                #    input_temp.append(os.path.join(outputdir, app_output_dir, app_outfile ))
                #    input_pdbs = sorted(input_temp)

    return_output_dict = {}

    MutagenesisID_list, mutations, PredID_Details, df, REU_list = Fetch_Mutant_ID(predID)

    for Wildtype_PDB_ID, Mutant_PDB_ID, MutagenesisID_index in MutagenesisID_list:
        fresh_wt_pdb, tmp_wt_pdb, fresh_mut_pdb, tmp_mut_pdb, residue_maps, wt_to_mut_chains = Generate_PDBs_and_Resmaps(PredID_Details, Mutant_PDB_ID, MutagenesisID_index, df)

        # PCA ANAYLSIS STUFF
        # PCA_Analysis(predID, input_pdbs, tmp_mut_pdb)

        # Get prody coordinates from RosettaOut Ensemble structures
        point_mutants, point_acceptable_atoms = mutant_coordinates(input_pdbs, mutations, residue_maps, wt_to_mut_chains, tmp_mut_pdb, tmp_wt_pdb, input_type = 'RosettaOut')
        neighbors = find_neighbors(mutations, fresh_wt_pdb, 8)
        neighborhood, neighborhood_acceptable_atoms = neighborhood_coordinates(neighbors, input_pdbs, residue_maps, wt_to_mut_chains, tmp_mut_pdb, tmp_wt_pdb, input_type='RosettaOut')
        global_ca, global_acceptable_atoms = global_ca_coordinates(input_pdbs, tmp_mut_pdb, tmp_wt_pdb, residue_maps, wt_to_mut_chains, input_type = 'RosettaOut')

        # Get prody coordinates from reference mutant crystal structure
        mut_pdb_global = global_ca_coordinates(input_pdbs, tmp_mut_pdb, tmp_wt_pdb, residue_maps, wt_to_mut_chains, input_type = 'Mutant PDB')
        mut_pdb_neighbors = neighborhood_coordinates(neighbors, input_pdbs, residue_maps, wt_to_mut_chains, tmp_mut_pdb, tmp_wt_pdb, input_type='Mutant PDB')
        mut_pdb_points = mutant_coordinates(input_pdbs, mutations, residue_maps, wt_to_mut_chains, tmp_mut_pdb, tmp_wt_pdb, input_type = 'Mutant PDB')

        # Get prody coordinates for Mutant and WT reference crystal structure backbone
        #### wt_global_bb
        wt_pdb = prody.parsePDB(tmp_wt_pdb)
        atom_list = []
        for atom in wt_pdb.select('name CA').getHierView().getAtoms():
            if atom.getChid() + str(atom.getResnum()) in global_acceptable_atoms['RosettaOut']:
                atom_list.append(str(atom.getIndex()))
        wt_global_bb = [wt_pdb.select('index ' + ' '.join(atom_list))]
        global_ca_bb = mut_pdb_global

        #### wt_neighbors_bb
        atom_list = []
        for atom in wt_pdb.select('name C CA N').getHierView().getAtoms():
            if (atom.getChid(), atom.getResnum(), atom.getName()) in [(waffle[0], waffle[2], waffle[3]) for waffle in neighborhood_acceptable_atoms]:
                atom_list.append(str(atom.getIndex()))
        wt_neighbors_bb = [wt_pdb.select('index ' + ' '.join(atom_list))]
        neighborhood_bb = [asdf.select('name C CA N') for asdf in mut_pdb_neighbors]

        #### wt_points_bb
        wt_points_bb = {}
        point_mutants_bb = {}
        for mutant_residue in point_mutants:
            wt_points_bb[mutant_residue] = [wt_pdb.select('name C CA N and chain %s and resnum %s' %(mutant_residue[:1], mutant_residue[1:]))]
            point_mutants_bb[mutant_residue] = [waffle.select('name C CA N') for waffle in mut_pdb_points[mutant_residue]]

        # Add stuff to Output Dictionary
        return_output_dict['%s : %s' %(Wildtype_PDB_ID, Mutant_PDB_ID)] = {}
        #### Calculate RosettaOut - Mutant Reference RMSDs
        return_output_dict['%s : %s' %(Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant - Point_Mutant RMSD'] = rmsd( pyrmsd_calc, point_mutants, mut_pdb_points)
        return_output_dict['%s : %s' %(Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant - Neighborhood RMSD'] = rmsd(pyrmsd_calc, neighborhood, mut_pdb_neighbors)
        return_output_dict['%s : %s' %(Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant - Global RMSD'] = rmsd(pyrmsd_calc, global_ca, mut_pdb_global)

        # #### Calculate RosettaOut - WT Reference RMSDs
        # return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['WT_BB - Point_Mutant RMSD'] = rmsd(pyrmsd_calc, point_mutants_bb, wt_points_bb)
        # return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['WT_BB - Neighborhood RMSD'] = rmsd(pyrmsd_calc, neighborhood_bb, wt_neighbors_bb)
        # return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['WT_BB - Global RMSD'] = rmsd(pyrmsd_calc, global_ca_bb, wt_global_bb)

        #### RosettaOut Rossetta Energy Scores
        return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant Complex REUs'] = REU_list

        #### Mutant PDB Point Mutant X Angles

        # DEBUGGING
        pprint.pprint(residue_maps)
        pprint.pprint(wt_to_mut_chains)

        mutations_mut_num = []
        for mutation in mutations:
            print mutation
            print wt_to_mut_chains[mutation[1]]
            print residue_maps[(mutation[1], wt_to_mut_chains[mutation[1]])]['%s %s ' % (mutation[1], ('   ' + str(mutation[0]))[-3:])].split()[1]
            mutations_mut_num.append([residue_maps[(mutation[1], wt_to_mut_chains[mutation[1]])]['%s %s ' % (mutation[1], ('   ' + str(mutation[0]))[-3:])].split()[1],
                                      wt_to_mut_chains[mutation[1]],
                                      mutation[2],
                                      mutation[3]]
                                     )

        return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant PDB X Angles'] = chi_angles([tmp_mut_pdb], mutations_mut_num)

        return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['WT-Mutant BackBone RMSDs'] = {}
        return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['WT-Mutant BackBone RMSDs']['Global'] = rmsd(pyrmsd_calc, global_ca_bb, wt_global_bb)
        return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['WT-Mutant BackBone RMSDs']['Neighborhood'] = rmsd(pyrmsd_calc, neighborhood_bb, wt_neighbors_bb)
        return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['WT-Mutant BackBone RMSDs']['Point_Mutant'] = rmsd(pyrmsd_calc, point_mutants_bb, wt_points_bb)

        chi_angles_output = chi_angles(input_pdbs, mutations)
        #### X angle output

        # print return_output_dict['%s : %s' %(Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant - Point_Mutant RMSD'][0]
        for point_mutant in return_output_dict['%s : %s' %(Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant - Point_Mutant RMSD']:
            if chi_angles_output != {}:
                return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant - Point_Mutant RMSD'][point_mutant]['X Angles'] = chi_angles_output[point_mutant]
            else:
                return_output_dict['%s : %s' % (Wildtype_PDB_ID, Mutant_PDB_ID)]['Mutant - Point_Mutant RMSD'][point_mutant]['X Angles'] = 'No X Angles!'

    return return_output_dict

def unzip_to_tmp_dir(prediction_id, zip_file):
    tmp_dir = tempfile.mkdtemp(prefix='unzip_to_tmp_')
    unzip_path = os.path.join(tmp_dir, '%d-ddg' % prediction_id)
    os.makedirs(unzip_path)
    with zipfile.ZipFile(zip_file, 'r') as job_zip:
        job_zip.extractall(unzip_path)
    return tmp_dir

def multiprocessing_stuff(predID):
    # Define things
    PredID_output_dict = {}

    try:
        my_tmp_dir = unzip_to_tmp_dir(predID, '%s.zip' % predID)
        outdir = os.path.join(my_tmp_dir, '%s-ddg' % predID)
        print '*** %s has been unzipped!!! ***' % predID
        print "\n*** Calculating RMSDs for %s ***\n" % outdir

        try:
            PredID_output_dict[predID] = do_math( outdir, predID)
            print 'Calculations complete!!!'
        except Exception as fuuuu:
            PredID_output_dict = fuuuu
            print 'Failure %s : %s' % (predID, fuuuu)
            import traceback
            traceback.print_exc()

        shutil.rmtree(my_tmp_dir)
        return PredID_output_dict

    except IOError as MIA:
        print 'IO Failure %s: %s' % (predID, MIA)
        PredID_output_dict = {predID: 'Failure : %s' % MIA}
        return PredID_output_dict

    except:
        print 'General Failure: %s' % predID
        PredID_output_dict = {predID: 'General Failure'}
        return PredID_output_dict

def main():
    os.chdir('/kortemmelab/shared/DDG/ppijobs')

    # Old code, can technically still use this but it's just easier to add PredIDs to PredID_list using range() for loop
    # Grabs Prediction IDs from .csv file
    # csv_info = pd.read_csv('Location of data.csv for job_name')
    # PredID_list = []
    # for index, row in csv_info.iterrows():
    #     PredID_list.append(int(row[0].split()[0]))

    # CHANGE FOR EACH RUN
    job_name = 'ddg_analysis_type_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-score_method_Rescore-Talaris2014'
    # PredID_list = [94009, 94011, 94012, 94075, 94205, 94213, 94230, 94231, 94268, 94269, 94270, 94271, 94272, 94314, 94315, 94316, 94317, 94318, 94319, 94367, 94531, 94533, 94534, 94535, 94536, 94537, 94538, 94539, 94540, 94541, 94550, 94571, 94574, 94578, 94581, 94953, 94956, 94964, 94981, 95074, 95079, 95113, 95118, 95127, 95131, 95163]

    # job_name = 'zemu-psbrub_1.6-pv-nt50000-bruball'
    # PredID_list = [89049, 89051, 89052, 89115, 89245, 89253, 89270, 89271, 89308, 89309, 89310, 89311, 89312, 89354, 89355, 89356, 89357, 89358, 89359, 89407, 89571, 89573, 89574, 89575, 89576, 89577, 89578, 89579, 89580, 89581, 89590, 89611, 89614, 89618, 89621, 89993, 89996, 90004, 90021, 90114, 90119, 90153, 90158, 90167, 90171, 90203]

    # DEBUGGING
    PredID_list = [int(sys.argv[1])] #94535 has an RMSD of 22A for some reason...

    pool = multiprocessing.Pool(25)
    allmyoutput = pool.map(multiprocessing_stuff, PredID_list, 1)
    pool.close()
    pool.join()

    print allmyoutput

    print 'Dumping information to pickle'
    import pickle
    with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/StructuralMetrics-%s.pickle' % job_name, 'wb') as outfile:
        output_dict = {}
        for resultdict in allmyoutput:
            output_dict[resultdict.keys()[0]] = resultdict
        pickle.dump(output_dict, outfile, 0)

    # with open('Structural_metrics.pickle', 'rb') as outfile:
    #     b = pickle.load(outfile)

main()

# main(): uses pool to pass individual Prediction IDs to multiprocessing_stuff()
# multiprocessing_stuff(): Unzips PredID zip into a temp dir, runs do_math() function, then deletes temp dir
