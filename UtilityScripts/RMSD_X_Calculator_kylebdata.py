# This script is intended to calculate the average pairwise RMSD between all Rosetta-generated mutant ensemble members
# as a measure of ensemble diversity

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
import pprint
import sys
import re
import pandas as pd
import os
import tempfile
import zipfile
import shutil
import json
import multiprocessing

from kddg.api import settings
sys_settings = settings.load()

#Return mutations for Prediction ID
def Fetch_PredID_Info(predID):
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

    pdb_filename = PredID_Details['Files']['Input'][0]['Filename']
    pdb_chains = re.sub('_|\.', ' ', pdb_filename).split()[1]
    raw_pdb = PredID_Details['Files']['Input'][0]['Content']
    fresh_pdb = strip_pdbs(raw_pdb, pdb_chains)

    return mutations, fresh_pdb

def strip_pdbs(raw_pdb, pdb_chains):
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w+b', delete=False) as temp_PDBFile:
        delete_me_later = temp_PDBFile.name
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

    with open(delete_me_later, mode='rb+') as temp_PDBFile:
        parser = PDBParser(PERMISSIVE=1)
        fresh_pdb = parser.get_structure('Open', temp_PDBFile.name)

    os.remove(delete_me_later)

    return fresh_pdb

#From Kyle's Finalize.py
def find_neighbors( pdb_path, predID, mutations, open_strct, neighbor_distance = 8.0):
    #mutations = mutations_resfile(filenum_dir)
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

#Kyle's RMSDWrapper.py (Spread out all over the place now)
def rmsd( pyrmsd_calc, coordinates ):
    def generate_out_dict(rmsd):
        from pyRMSD.condensedMatrix import CondensedMatrix
        out_dict_temp = {}
        out_dict_temp['Mean'] = CondensedMatrix(rmsd).calculateMean()
        out_dict_temp['Variance'] = CondensedMatrix(rmsd).calculateVariance()
        out_dict_temp['Skewness'] = CondensedMatrix(rmsd).calculateSkewness()
        out_dict_temp['Kurtosis'] = CondensedMatrix(rmsd).calculateKurtosis()
        return out_dict_temp

    if type(coordinates) is dict:
        rmsd_list = []
        for mutres, coord_set in coordinates.iteritems():
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(pyrmsd_calc, coord_set)
            rmsd = calculator.pairwiseRMSDMatrix()
            out_dict = generate_out_dict(rmsd)
            rmsd_list.append([mutres, out_dict])
        return rmsd_list
    else:
        # rmsd_matrix = pyRMSD.matrixHandler.MatrixHandler().createMatrix(coordinates, pyrmsd_calc)
        # print rmsd_matrix
        # print type(coordinates)
        # avg_rmsd = average_rmsd_calculator(rmsd_matrix)

        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(pyrmsd_calc, coordinates)
        rmsd = calculator.pairwiseRMSDMatrix()
        out_dict = generate_out_dict(rmsd)
        return out_dict

def average_rmsd_calculator(rmsd_matrix):
    avg_rmsd_list = rmsd_matrix.get_data()
    avg_rmsd = sum(avg_rmsd_list) / len(avg_rmsd_list)
    return avg_rmsd

def global_ca_coordinates(input_pdbs):
    coordinates = np.array( [prody.parsePDB(input_pdb).select('calpha').getCoords() for input_pdb in input_pdbs] )
    return coordinates    
    
def neighborhood_coordinates(neighbors, reference, input_pdbs):
    temp = []

    for input_pdb in input_pdbs:
        if 'WT.' not in input_pdb:
            atom_list = []
            hv = prody.parsePDB(input_pdb).getHierView()
            res_list = [hv[neighbor[1], neighbor[0][1]] for neighbor in neighbors]
            for res in res_list:
                for atom in res:
                    if atom.getElement() != 'H':
                        atom_list.append(atom.getCoords())
                    temp.append(atom_list)
        
    coordinates = np.asarray(temp)
    return coordinates
        
def mutant_coordinates(input_pdbs, predID, mutations, PDBFile_from_database):
    mutation_dict = {}

    for mutation in mutations:
        temp = []
        temp_nparray = []
        
        for input_pdb in input_pdbs: 
            if 'WT.' not in input_pdb:
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
    return mutation_dict

def bin_me(templist):
    from collections import Counter
    bins= {}
    floor_list = [X // 30 for X in templist]
    X_counts = Counter(floor_list)
    bin_ranges = {int(-6.0) : '-180 <= x < -150',
                  int(-5.0) : '-150 <= x < -120',
                  int(-4.0) : '-120 <= x < -90',
                  int(-3.0) : '-90 <= x < -60',
                  int(-2.0) : '-60 <= x <-30',
                  int(-1.0) : '-30 <= x < 0',
                  int(0) : '0 <= x < 30',
                  int(1.0) : '30 <= x < 60',
                  int(2.0) : '60 <= x < 90',
                  int(3.0) : '90 <= x < 120',
                  int(4.0) : '120 <= x < 150',
                  int(5.0) : '150 <= x <180',
                  int(6.0) : 'x = 180'}
    readable_X_counts = dict((bin_ranges[key], value) for (key, value) in X_counts.items())
    return readable_X_counts
            
def chi_angles(input_pdbs, predID, mutations, PDBFile_from_database):
    
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
        if mutation[3] == 'A' or mutation[3] == 'G':
            xangles_dict['%s%s' %(mutation[0], mutation[1])] = "Mutation is %s!" %mutation[3]
        else:
            xangles_dict['%s%s' %(mutation[0], mutation[1])] = {}
            x1_templist = []
            x2_templist = []
            for input_pdb in input_pdbs:
                if 'WT.pdb_' not in input_pdb:
                    p = prody.parsePDB(input_pdb)
                    hv = p.getHierView()

                    current_res = hv[mutation[1], int(mutation[0])]
                    
                    if current_res.getResname() in chi1_dict.keys():
                        current_res.select('name %s' %chi1_dict[current_res.getResname()])
                        chi1 = prody.measure.measure.calcDihedral(current_res.select('name %s' %chi1_dict[current_res.getResname()][0]),
                                                              current_res.select('name %s' %chi1_dict[current_res.getResname()][1]),
                                                              current_res.select('name %s' %chi1_dict[current_res.getResname()][2]),
                                                              current_res.select('name %s' %chi1_dict[current_res.getResname()][3]))
                        x1_templist.append(chi1[0])
                        
                    if current_res.getResname() in chi2_dict.keys():
                        chi2 = prody.measure.measure.calcDihedral(current_res.select('name %s' %chi2_dict[current_res.getResname()][0]),
                                                              current_res.select('name %s' %chi2_dict[current_res.getResname()][1]),
                                                              current_res.select('name %s' %chi2_dict[current_res.getResname()][2]),
                                                              current_res.select('name %s' %chi2_dict[current_res.getResname()][3]))
                        x2_templist.append(chi2[0])

            x1_bins = bin_me(x1_templist)
            xangles_dict['%s%s' % (mutation[0], mutation[1])]['X1'] = x1_bins

            if x2_templist != []:
                x2_bins = bin_me(x2_templist)
                xangles_dict['%s%s' % (mutation[0], mutation[1])]['X2'] = x2_bins

            #Outputs Ramachandran-like plots for X1 and X2, Violin plots (in progress) if only X1 is present
            #import pandas as pd
            #import seaborn as sns
            #import matplotlib.pyplot as plt
            
            #x1 = pd.Series(x1_templist, name = 'X1')
            #if x2_templist != []:
            #    x2 = pd.Series(x2_templist, name = 'X2')
            #    mrplotty = sns.regplot(x1, x2, fit_reg = False)
            #    plt.savefig('%s%splot.pdf'  %(mutation[0], mutation[1]))

    return xangles_dict

#Action!!!
def do_math( reference, outputdir, predID):
    # Use CUDA for GPU calculations, if avialable
    if 'QCP_CUDA_MEM_CALCULATOR' in availableCalculators():
        pyrmsd_calc = 'QCP_CUDA_MEM_CALCULATOR'
    else:
        pyrmsd_calc = 'QCP_SERIAL_CALCULATOR'
        
    input_temp = []

    for app_output_dir in os.listdir(outputdir):
        for app_outfile in os.listdir(os.path.join(outputdir, app_output_dir)):
            #Mutant Structures
            if app_outfile == 'min_cst.mutate_bkrb_min_cst_0001_0001.pdb.gz':
                input_temp.append(os.path.join(outputdir, app_output_dir, app_outfile))
                input_pdbs = sorted(input_temp)
            # WT Structures
            #if app_outfile == 'min_cst.repack-wt_bkrb_min_cst_0001_0001.pdb.gz':
            #    input_temp.append(os.path.join(outputdir, app_output_dir, app_outfile ))
            #    input_pdbs = sorted(input_temp)

    mutations, fresh_pdb = Fetch_PredID_Info(predID)

    point_mutants = mutant_coordinates(input_pdbs, predID, mutations, fresh_pdb)
    neighbors = find_neighbors( reference, predID, mutations, fresh_pdb, 8)
    neighborhood = neighborhood_coordinates(neighbors, reference, input_pdbs)
    global_ca = global_ca_coordinates(input_pdbs)

    return_output_dict = {}

    return_output_dict['Point Mutant RMSDs'] = rmsd( pyrmsd_calc, point_mutants)
    return_output_dict['Neighborhood RMSD'] = rmsd(pyrmsd_calc, neighborhood)
    return_output_dict['Global RMSD'] = rmsd(pyrmsd_calc, global_ca)
    chi_angles_output = chi_angles(input_pdbs, predID, mutations, fresh_pdb)
    if chi_angles_output != {}:
        return_output_dict['X angles'] = chi_angles_output
    else:
        print 'No X angles!'
    return return_output_dict

def unzip_to_tmp_dir(prediction_id, zip_file):
    tmp_dir = tempfile.mkdtemp(prefix='unzip_to_tmp_')
    unzip_path = os.path.join(tmp_dir, '%d-ddg' % prediction_id)
    os.makedirs(unzip_path)
    with zipfile.ZipFile(zip_file, 'r') as job_zip:
        job_zip.extractall(unzip_path)
    return tmp_dir

def multiprocessing_stuff(predID):
    def dict_to_json(PredID_output_dict, predID):
        with open('/kortemmelab/home/james.lucas/Structural_metrics_Outfiles/%s.txt' % predID, 'a') as outfile:
            json.dump(PredID_output_dict, outfile)
    try:
        my_tmp_dir = unzip_to_tmp_dir(predID, '%s.zip' %predID)

        print '*** %s has been unzipped!!!***' %predID

        # Define things
        PredID_output_dict = {}
        outdir = os.path.join(my_tmp_dir, '%s-ddg' % predID)

        print "\n***Calculating RMSDs for %s***\n" % outdir
        reference = '../data/%s/%s'  # %(outdir, mypdb_ref)

        try:
            PredID_output_dict[predID] = do_math(reference, outdir, predID)
            print 'Calculations complete!!!'
        except Exception as fuuuu:
            PredID_output_dict[predID] = fuuuu
            print 'Failure %s : %s' % (predID, fuuuu)
            import traceback
            traceback.print_exc()

        shutil.rmtree(my_tmp_dir)
        dict_to_json(PredID_output_dict, predID)
        return PredID_output_dict

    except IOError as MIA:
        print 'IO Failure %s: %s' %(predID, MIA)
        PredID_output_dict = {'%s' %predID: '%s' %MIA}
        dict_to_json(PredID_output_dict, predID)
        return PredID_output_dict

    except:
        print 'General Failure: %s' %predID
        PredID_output_dict = {'%s' %predID : 'General Failure'}
        dict_to_json(PredID_output_dict, predID)
        import traceback
        traceback.print_exc()
        return PredID_output_dict

def main():
    os.chdir('/kortemmelab/shared/DDG/ppijobs')
    csv_info = pd.read_csv('/kortemmelab/home/james.lucas/zemu-psbrub_1.6-pv-1000-lbfgs_IDs.csv')

    # Grab Prediction IDs from .csv file
    # PredID_list = []
    # for index, row in csv_info.iterrows():
    #     PredID_list.append(int(row[0].split()[0]))

    PredID_list = []
    # for asdf in range(67088, 68327): #zemu-psbrub_1.6-pv-1000-lbfgs
    for asdf in range(94009,95248): #ddg_analysis_type_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-score_method_Rescore-Talaris2014
        PredID_list.append(asdf)

    #DEBUGGING
    #PredID_list = [67619, 67611, 67614]
    PredID_list = [94009]

    pool = multiprocessing.Pool(25)
    allmyoutput = pool.map( multiprocessing_stuff, PredID_list, 1)
    pool.close()
    pool.join()

    print allmyoutput

    # print 'Dumping information to Structural_metrics.txt'
    #
    # with open('/kortemmelab/home/james.lucas/Structural_metrics.txt', 'a') as outfile:
    #     for resultdict in allmyoutput:
    #         json.dump(resultdict, outfile)

main()
