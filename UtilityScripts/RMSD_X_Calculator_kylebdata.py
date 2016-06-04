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

#Return mutations for Prediction ID
def PredID_to_mutations(predID):
    # List of dictionaries for Resfile mutation info
    sys.path.insert(0,
                    '/kortemmelab/home/james.lucas')  # this should point to the directory above your ddg repo checkout
    from ddg.ddglib.ppi_api import get_interface_with_config_file as get_ppi_interface

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
def find_neighbors(filenum_dir, pdb_path, predID, neighbor_distance = 8.0):
    #mutations = mutations_resfile(filenum_dir)
    mutations, open_strct = PredID_to_mutations(predID)

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
def rmsd(input_pdbs, pyrmsd_calc, coordinates):
    if type(coordinates) is dict:
        rmsd_list = []
        for mutres, coord_set in coordinates.iteritems():
            rmsd_matrix = pyRMSD.matrixHandler.MatrixHandler().createMatrix(coord_set, pyrmsd_calc)
            avg_rmsd = average_rmsd_calculator(rmsd_matrix)
            rmsd_list.append([mutres, avg_rmsd])
        return rmsd_list
    else:
        rmsd_matrix = pyRMSD.matrixHandler.MatrixHandler().createMatrix(coordinates, pyrmsd_calc)
        avg_rmsd = average_rmsd_calculator(rmsd_matrix)
        return avg_rmsd
    
def average_rmsd_calculator(rmsd_matrix):
    rmsd_sqform = scipy.spatial.distance.squareform( rmsd_matrix.get_data() )
    avg_rmsd_list = []
    for rmsd_row in enumerate(rmsd_sqform):
        for rmsd_value in enumerate(rmsd_row[1]):
            if rmsd_value[0] > rmsd_row[0]:
                avg_rmsd_list.append(rmsd_value[1])
    avg_rmsd = sum(avg_rmsd_list)/len(avg_rmsd_list)
    return avg_rmsd
    
def global_ca_rms(input_pdbs):
    coordinates = np.array( [prody.parsePDB(input_pdb).select('calpha').getCoords() for input_pdb in input_pdbs] )
    return coordinates    
    
def neighborhood_rms(neighbors, reference, input_pdbs):
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
                    else:
                        pass
                    temp.append(atom_list)
        
    coordinates = np.asarray(temp)
    return coordinates
        
#Calculates sidechain rmsd        
def mutant_rms(datadir, input_pdbs, predID):
    #mutations = read_mutations_resfile(datadir)
    mutations, PDBFile_from_database = PredID_to_mutations(predID)
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

#Bins X-angles
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
            
def chi_angles(datadir, input_pdbs, predID):
    mutations, PDBFile_from_database = PredID_to_mutations(predID)
    
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

            #Bin Xangles and add to bin dicts
            x1_bins = bin_me(x1_templist)
            if x2_templist != []:
                x2_bins = bin_me(x2_templist)

            #Add things to xangles_dict!
            xangles_dict['%s%s' %(mutation[0], mutation[1])]['X1'] = x1_bins
            xangles_dict['%s%s' %(mutation[0], mutation[1])]['X2'] = x2_bins

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
def do_math(datadir, reference, outputdir, predID):
    # Use CUDA for GPU calculations, if avialable
    if 'QCP_CUDA_MEM_CALCULATOR' in availableCalculators():
        pyrmsd_calc = 'QCP_CUDA_MEM_CALCULATOR'
    else:
        pyrmsd_calc = 'QCP_SERIAL_CALCULATOR'
        
    input_temp = []

    for app_output_dir in os.listdir(outputdir):
        for app_outfile in os.listdir(app_output_dir):
            #Mutant Structures
            if app_outfile == 'min_cst.mutate_bkrb_min_cst_0001_0001.pdb.gz':
                input_temp.append(os.path.join(outputdir, app_output_dir, app_outfile))
                input_pdbs = sorted(input_temp)
            # WT Structures
            #if app_outfile == 'min_cst.repack-wt_bkrb_min_cst_0001_0001.pdb.gz':
            #    input_temp.append(os.path.join(outputdir, app_output_dir, app_outfile ))
            #    input_pdbs = sorted(input_temp)

    neighbors = find_neighbors(datadir, reference, predID, 8)
    point_mutants = mutant_rms(datadir, input_pdbs, predID)
    neighborhood = neighborhood_rms(neighbors, reference, input_pdbs)
    global_ca = global_ca_rms(input_pdbs)

    return_output_dict = {}

    return_output_dict['Point Mutant RMSDs'] = rmsd(input_pdbs, pyrmsd_calc, point_mutants)
    return_output_dict['Neighborhood RMSD'] = rmsd(input_pdbs, pyrmsd_calc, neighborhood)
    return_output_dict['Global RMSD'] = rmsd(input_pdbs, pyrmsd_calc, global_ca)
    chi_angles_output = chi_angles(datadir, input_pdbs, predID)
    if chi_angles_output != {}:
        return_output_dict['X angles'] = chi_angles_output
    else:
        print 'No X angles!'
    return return_output_dict

def main(predID, my_tmp_dir):
    #Define things
    output_dict = {}
    outdir = os.path.join(my_tmp_dir, '%s-ddg' %predID)
    datadir = '../data/%s' #%outdir
    #Change to root of output directory
    #os.chdir('/kortemmelab/home/james.lucas/160412-kyleb_jl-brub-rscr-v2/DDG_Zemu_v2_output-Sum_DDG_only')
    os.chdir(outdir)
    
    print "\n***Calculating RMSDs for %s***\n" %outdir

    reference = '../data/%s/%s' #%(outdir, mypdb_ref)

    # #DEBUGGING ONLY
    # PredID_to_mutations(predID)

    output_dict[predID] = do_math(datadir, reference, outdir, predID)

    return output_dict

def unzip_to_tmp_dir(prediction_id, zip_file):
    tmp_dir = tempfile.mkdtemp(prefix='unzip_to_tmp_')
    unzip_path = os.path.join(tmp_dir, '%d-ddg' % prediction_id)
    os.makedirs(unzip_path)
    with zipfile.ZipFile(zip_file, 'r') as job_zip:
        job_zip.extractall(unzip_path)
    return tmp_dir

def multiprocessing_stuff(predID):
    my_tmp_dir = unzip_to_tmp_dir(predID, '%s.zip' %predID)
    print '*** %s has been unzupped!!!***' %predID
    try:
        PredID_output_dict = main(predID, my_tmp_dir)
    except:
        PredID_output_dict = {predID}
        print 'ERROR: %s' %predID
    shutil.rmtree(my_tmp_dir)
    return PredID_output_dict

def asdfasdf():
    os.chdir('/kortemmelab/shared/DDG/ppijobs')
    csv_info = pd.read_csv('/kortemmelab/home/james.lucas/zemu-psbrub_1.6-pv-1000-lbfgs_IDs.csv')

    PredID_list = []
    for index, row in csv_info.iterrows():
        PredID_list.append(int(row[0].split()[0]))

    #DEBUGGING

    pool = multiprocessing.Pool(2)
    allmyoutput = pool.map( multiprocessing_stuff, PredID_list, 1)
    pool.close()
    pool.join()


    print allmyoutput

    with open('/kortemmelab/home/james.lucas/Structural_metrics.txt', 'a') as outfile:
        for resultdict in allmyoutput:
            json.dump(resultdict, outfile)

asdfasdf()
