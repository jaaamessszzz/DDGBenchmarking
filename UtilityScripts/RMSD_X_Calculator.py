from Bio.PDB import *
import Bio.PDB
import os
import pyRMSD
from pyRMSD.matrixHandler import MatrixHandler
import pyRMSD.RMSDCalculator
from pyRMSD.availableCalculators import availableCalculators
import prody
import numpy as np
import scipy.spatial.distance

#From Kyle's Finalize.py
def read_mutations_resfile(filenum_dir):
    resfile = os.path.join(filenum_dir, 'mutations_repack.resfile')
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

#From Kyle's Finalize.py
def find_neighbors(filenum_dir, pdb_path, neighbor_distance = 8.0):
    mutations = read_mutations_resfile(filenum_dir)
    parser = PDBParser(PERMISSIVE=1)
    open_strct = parser.get_structure('Open', pdb_path)

    # There should only be one model in PDB file
    num_models = 0
    for model in open_strct.get_models():
        num_models += 1
    assert( num_models == 1 )

    chain_list = [chain.get_id() for chain in open_strct[0].get_chains()]
    neighbors = set()
    for mutation in mutations:
        res_id, chain_id, pikaa, mut_aa = mutation
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
                    try:
                        dist = mut_residue['CA'] - open_strct[0][chain][residue]['CA']
                        if dist < neighbor_distance:
                            neighbors.add( (residue, chain) )
                    except KeyError:
                        pass

    return neighbors

#Kyle's RMSDWrapper.py (Spread out al over the place now)
def rmsd(input_pdbs, pyrmsd_calc, coordinates):
    if type(coordinates) is dict:
        rmsd_list = []
        for mutres, coord_set in coordinates.iteritems():
            rmsd_matrix = pyRMSD.matrixHandler.MatrixHandler().createMatrix(coord_set, pyrmsd_calc)
            rmsd_list.append([mutres, scipy.spatial.distance.squareform( rmsd_matrix.get_data() )])
        return rmsd_list
    else:
        rmsd_matrix = pyRMSD.matrixHandler.MatrixHandler().createMatrix(coordinates, pyrmsd_calc)
        return [scipy.spatial.distance.squareform( rmsd_matrix.get_data() )]
    
def global_ca_rms(input_pdbs):
    coordinates = np.array( [prody.parsePDB(input_pdb).select('calpha').getCoords() for input_pdb in input_pdbs] )
    return coordinates    
    
def neighborhood_rms(neighbors, reference, input_pdbs):
    temp = []

    for input_pdb in input_pdbs:
        atom_list = []
        hv = prody.parsePDB(input_pdb).getHierView()
        res_list = sorted( [hv[neighbor[1], neighbor[0][1]] for neighbor in neighbors] )
        for res in res_list:
            for atom in res:
                if atom.getElement() != 'H':
                    print atom
                    atom_list.append(atom.getCoords())
                else:
                    pass
        temp.append(atom_list)
        
    coordinates = np.asarray(temp)
    return coordinates
        
#Calculates sidechain rmsd        
def mutant_rms(datadir, input_pdbs):
    mutations = read_mutations_resfile(datadir)
    mutation_dict = {}
    
    for mutation in mutations:
        temp = []
        temp_nparray = []
        
        for input_pdb in input_pdbs: 
            atom_list = []
            hv = prody.parsePDB(input_pdb).getHierView()
            res_list = [hv[mutation[1], int(mutation[0])]]
            for res in res_list:
                for atom in res:
                    if atom.getElement() != 'H':
                        atom_list.append(atom.getCoords())
                    else:
                        pass
            temp.append(atom_list)
        
        temp_nparray = np.asarray(temp)
        
        mutation_dict['%s%s' %(mutation[1], mutation[0])] = temp_nparray
    print mutation_dict.keys()
    return mutation_dict

def chi_angles(datadir, input_pdbs):
    mutations = read_mutations_resfile(datadir)
    
    #http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
    chi1_dict = {'ARG':['N','CA','CB','CG'],
                 'ASN':['N','CA','CB','CG'],
                 'ASP':['N','CA','CB','CG'],
                 'CYS':['N','CA','CB','SG'],
                 'GLN':['N','CA','CB','CG'],
                 'GLU':['N','CA','CB','CG'],
                 'HIS':['N','CA','CB','CG'],
                 'ILE':['N','CA','CB','CG1'],
                 'LEU':['N','CA','CB','CG'],
                 'LYS':['N','CA','CB','CG'],
                 'MET':['N','CA','CB','CG'],
                 'PHE':['N','CA','CB','CG'],
                 'PRO':['N','CA','CB','CG'],
                 'SER':['N','CA','CB','OG'],
                 'THR':['N','CA','CB','OG1'],
                 'TRP':['N','CA','CB','CG'],
                 'TYR':['N','CA','CB','CG'],
                 'VAL':['N','CA','CB','CG1']}
    chi2_dict = {'ARG':['CA','CB','CG','CD'],
                 'ASN':['CA','CB','CG','OD1'],
                 'ASP':['CA','CB','CG','OD1'],
                 'GLN':['CA','CB','CG','CD'],
                 'GLU':['CA','CB','CG','CD'],
                 'HIS':['CA','CB','CG','ND1'],
                 'ILE':['CA','CB','CG1','CD'],
                 'LEU':['CA','CB','CG','CD1'],
                 'LYS':['CA','CB','CG','CD'],
                 'MET':['CA','CB','CG','SD'],
                 'PHE':['CA','CB','CG','CD1'],
                 'PRO':['CA','CB','CG','CD'],
                 'TRP':['CA','CB','CG','CD1'],
                 'TYR':['CA','CB','CG','CD1']}
    
    xangles_dict = {}
    
    for mutation in mutations:
        if mutation[3] == 'A' or 'G':
            continue
        else:
            xangles_dict['%s%s' %(mutation[0], mutation[1])] = {}
            counter = 1
            for input_pdb in input_pdbs:
                p = prody.parsePDB(input_pdb)
                hv = p.getHierView()

                current_res = hv[mutation[1], int(mutation[0])]
                templist = []

                if current_res.getResname() == 'ALA' or 'GLY':
                    continue

                else:
                    if current_res.getResname() in chi1_dict.keys():
                        current_res.select('name %s' %chi1_dict[current_res.getResname()])
                        chi1 = calcDihedral(current_res.select('name %s' %chi1_dict[current_res.getResname()][0]),
                                            current_res.select('name %s' %chi1_dict[current_res.getResname()][1]),
                                            current_res.select('name %s' %chi1_dict[current_res.getResname()][2]),
                                            current_res.select('name %s' %chi1_dict[current_res.getResname()][3]))
                        templist.append(chi1[0])

                    if current_res.getResname() in chi2_dict.keys():
                        chi2 = calcDihedral(current_res.select('name %s' %chi2_dict[current_res.getResname()][0]),
                                            current_res.select('name %s' %chi2_dict[current_res.getResname()][1]),
                                            current_res.select('name %s' %chi2_dict[current_res.getResname()][2]),
                                            current_res.select('name %s' %chi2_dict[current_res.getResname()][3]))
                        templist.append(chi2[0])

                    xangles_dict['%s%s' %(mutation[0], mutation[1])]['%s' %counter] = sorted(templist)

                counter = counter + 1

            print xangles_dict
    return xangles_dict

#Action!!!
def main():
    os.chdir('/home/james.lucas/Rotation/DDGBenchmarks_Test/')
    
    #Define things
    datadir = 'TestJobs/data/59648/'
    reference = 'TestJobs/data/59648/1TM1_EI.pdb'
    outputdir = 'TestJobs/output/59648/'
    
    # Use CUDA for GPU calculations, if avialable
    if 'QCP_CUDA_MEM_CALCULATOR' in availableCalculators():
        pyrmsd_calc = 'QCP_CUDA_MEM_CALCULATOR'
    else:
        pyrmsd_calc = 'QCP_SERIAL_CALCULATOR'
    
    input_temp = []
    for i in os.listdir(outputdir):
        if i.endswith('.pdb'):
            input_temp.append(os.path.join(outputdir, i ))
    input_pdbs = sorted(input_temp)
    
    neighbors = find_neighbors(datadir, reference, 8)
    
    point_mutants = mutant_rms(datadir, input_pdbs)
    #neighborhood = neighborhood_rms(neighbors, reference, input_pdbs)
    #global_ca = global_ca_rms(input_pdbs)
    
    asdf = rmsd(input_pdbs, pyrmsd_calc, point_mutants)
    
    for i in asdf:
        print i
        
main()