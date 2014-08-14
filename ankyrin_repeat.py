import sys
import os
import datetime
import string
import random
import glob
import shutil
import json
import itertools
sys.path.insert(0, "..")

import numpy as np
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt

from tools import colortext
from tools.bio.pdb import PDB
from tools.bio.basics import Mutation, ChainMutation, generate_all_combinations_of_mutations
from tools.bio.alignment import ScaffoldModelChainMapper
from tools.fs.fsio import read_file, write_file
from tools.process import Popen
from ddglib import ddgdbapi, dbapi

DDGdb = ddgdbapi.ddGDatabase()
ddG_connection = dbapi.ddG()

all_wildtype_mutations = '''
# All mutations are on chain A
A83P, I86L, S87A
L95I, I98V
A83P, I86L, S87A, L95I, I98V
'''

def add_pdb_file(filepath, pdb_id):
    existing_pdb = DDGdb.execute_select('SELECT ID FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))

    if not existing_pdb:
        pdb_contents = read_file(filepath)
        p = PDB(pdb_contents)

        fasta = []
        for c, sequence in p.atom_sequences.iteritems():
            fasta.append('>%s:%s|PDBID|CHAIN|SEQUENCE' % (pdb_id.replace(':', '_'), c))
            fasta.append(str(sequence))
        fasta = '\n'.join(fasta)

        d = {
            'ID' : pdb_id,
            'FileSource' : 'Biosensor project',
            'Content' : read_file(filepath),
            'FASTA' : fasta,
            'Resolution' : None,
            'Techniques' : 'Rosetta model',
            'BFactors' : '',
            'Publication' : None
        }
        DDGdb.insertDictIfNew('PDBFile', d, ['ID'])

def add_mutant(pdb_ID, mutant_mutations):
    '''Use this function to add one set of mutations ON THE SAME CHAIN (i.e. corresponding to one mutant) to the database.
       todo: generalize this to allow different chains
    '''
    chains = set([m.Chain for m in mutant_mutations])
    assert(len(chains) == 1)
    colortext.warning("Adding mutation: %s." % ', '.join(map(str, mutant_mutations)))
    ddG_connection.createDummyExperiment_ankyrin_repeat(pdb_ID, mutant_mutations, chains.pop())


def add_cluster_jobs(pdb_ID, PredictionSet, ProtocolID):
    colortext.printf("\nAdding any mutations for this structure which have not been queued/run in the %s prediction set." % PredictionSet, "lightgreen")

    KeepHETATMLines = False

    d = {
        'ID' : PredictionSet,
        'Status' : 'halted',
        'Priority' : 9,
        'BatchSize' : 40,
        'EntryDate' : datetime.datetime.now(),
    }
    DDGdb.insertDictIfNew('PredictionSet', d, ['ID'])

    # Determine the set of experiments to add
    ExperimentIDs = set([r['ID'] for r in DDGdb.execute_select('SELECT ID FROM Experiment WHERE PDBFileID=%s', parameters=(pdb_ID,))])
    ExperimentIDsInPredictionSet = set([r['ExperimentID'] for r in DDGdb.execute_select('SELECT ExperimentID FROM Prediction WHERE PredictionSet=%s', parameters=(PredictionSet,))])
    experiment_IDs_to_add = sorted(ExperimentIDs.difference(ExperimentIDsInPredictionSet))

    if experiment_IDs_to_add:
        colortext.printf("\nAdding %d jobs to the cluster queue." % len(experiment_IDs_to_add), "lightgreen")

        for experiment_ID in experiment_IDs_to_add:
            colortext.write('.', "lightgreen")
            ddG_connection.addPrediction(experiment_ID, None, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True)
    else:
        colortext.printf("\nAll jobs are already in the queue or have been run.", "lightgreen")
    print('')

def get_results(PredictionSet):
    return DDGdb.execute_select('''
SELECT Prediction.ID AS PredictionID, Prediction.ExperimentID, Experiment.PDBFileID, ExperimentMutations.FlattenedMutations, Prediction.Scores, TIMEDIFF(Prediction.EndDate, Prediction.StartDate) AS TimeTaken FROM Prediction INNER JOIN
(
  SELECT ExperimentID, GROUP_CONCAT(Mutation SEPARATOR ', ') AS FlattenedMutations FROM
  (
    SELECT ExperimentID, CONCAT(Chain, ' ', WildTypeAA, ResidueID, MutantAA) As Mutation FROM ExperimentMutation
  ) AS FlattenedMutation
  GROUP BY ExperimentID
) AS ExperimentMutations
ON Prediction.ExperimentID=ExperimentMutations.ExperimentID
INNER JOIN Experiment ON Prediction.ExperimentID=Experiment.ID
WHERE Prediction.PredictionSet=%s
ORDER BY Prediction.ExperimentID''', parameters=(PredictionSet,))

def analyze_results(PredictionSet, graph_filename, scoring_method, scoring_type):
    results = get_results(PredictionSet)
    sorted_results = {}
    for r in results:
        sorted_results[(json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], r['ExperimentID'])] = r
    count = 0

    set_of_mutations = set()

    for k, r in sorted(sorted_results.iteritems()):
        #if r['FlattenedMutations'].find('A E141L') != -1 and r['FlattenedMutations'].find('A S142A') != -1 and r['FlattenedMutations'].find('A L78Y') != -1:
        #    print('%f, %s' % (k[0], r['FlattenedMutations']))
        #if r['FlattenedMutations'].find('A W103M') != -1 and r['FlattenedMutations'].find('A F70Y') != -1:
        #    if r['FlattenedMutations'].find('A E141L') == -1 and r['FlattenedMutations'].find('A S142A') == -1 and r['FlattenedMutations'].find('A L78Y') == -1:
        #        print('%f, %s' % (k[0], r['FlattenedMutations']))

        if r['FlattenedMutations'].find('A W103M') != -1 and r['FlattenedMutations'].find('A F70Y') != -1:
            if r['FlattenedMutations'].find('A E141L') == -1 and r['FlattenedMutations'].find('A S142A') == -1 and r['FlattenedMutations'].find('A L78Y') == -1:
                #print('%f, %s' % (k[0], r['FlattenedMutations']))
                count += 1
        #A E141L, A S142A

        mutations = [m for m in map(string.strip, r['FlattenedMutations'].split(',')) if m]
        for m in mutations:
            set_of_mutations.add((int(m.split()[1][1:-1]), m))
        if r['FlattenedMutations'].find('A L78Y') == -1:
            print('%f, %s' % (k[0], r['FlattenedMutations']))
        #    #count += 1

    data = []
    pruned_data = []
    for k, r in sorted(sorted_results.iteritems()):
        line = []
        for m in sorted(set_of_mutations):
            if r['FlattenedMutations'].find(m[1]) != -1:
                line.append(1)
            else:
                line.append(0)
        data.append((json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], line))
        if r['FlattenedMutations'].find('A L78Y') == -1:
            pruned_data.append((json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], line))

    labels = [m[1] for m in sorted(set_of_mutations)]

    #hinton('all_mutants.png', labels, data)
    create_graph(graph_filename, labels, pruned_data, scoring_method, scoring_type)

def create_graph(filename, labels, data, scoring_method, scoring_type):
    matplotlib.rc('figure', figsize=(8.27, 20.69))

    x_values = []
    y_values = []
    ddg_values = []
    y = 0
    for line in data:
        x = 0
        y += 7
        w = line[0]
        plt.text(30, y, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=9)
        for point in line[1]:
            x += 3
            if point == 1:
                x_values.append(x)
                y_values.append(y)
                ddg_values.append(line[0])

    if len(data) > 140:
        plt.scatter(x_values, y_values, c=ddg_values, s=10, cmap=matplotlib.cm.jet, edgecolors='none', zorder=99)
    else:
        plt.scatter(x_values, y_values, c=ddg_values, s=50, cmap=matplotlib.cm.jet, edgecolors='none', zorder=99)

    plt.tight_layout(pad=2.08)
    plt.axis((0, 27, -5, (7 * len(data)) + 15))

    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        left='off',      # ticks along the left edge are off
        labelleft='off', # labels along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') # labels along the bottom edge are off

    x = 1.9
    for l in labels:
        plt.text(x, -12, l.split()[1], fontdict=None, withdash=True, fontsize=9)
        x += 3

    added_zero_line = False
    last_y_value = 0
    y = 0
    for line in data:
        x = 0
        y += 7
        plt.plot([1, 25], [y, y], color='#999999', linestyle='-', linewidth=0.1)
        if y % 21 == 7:
            if len(data) > 140:
                plt.text(25, y-3.5, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=6)
            else:
                plt.text(25, y-1.75, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=6)
        if not added_zero_line:
            if line[0] > 0:
                plt.plot([1, 25], [0.5 + ((y + last_y_value) / 2), 0.5 + ((y + last_y_value) / 2)], color='k', linestyle='-', linewidth=1)
                added_zero_line = True
            else:
                last_y_value = y

    if len(data) > 140:
        plt.text(25, y-3.5, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=6)
    else:
        plt.text(25, y-1.75, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=6)

    plt.colorbar()
    plt.title(r'$\Delta\Delta$G predictions for FPP biosensor mutants (%s.%s)' % (scoring_method.replace(',0A', '.0$\AA$').replace('_', ' '), scoring_type), fontdict=None)

    if len(data) > 140:
        plt.savefig(filename, dpi=800)
    else:
        plt.savefig(filename, dpi=400)


def extract_data(ddG_connection, output_dir, PredictionID):
    assert(os.path.exists(output_dir))
    colortext.message('here')
    archive = ddG_connection.getData(PredictionID)
    write_file(os.path.join(output_dir, '%d.zip' % PredictionID), archive, 'wb')
    colortext.message('Unzipping %d.zip' % PredictionID)
    p = Popen(output_dir, ['unzip', '%d.zip' % PredictionID])
    os.remove(os.path.join(output_dir, '%d.zip' % PredictionID))
    if p.errorcode != 0:
        raise colortext.Exception(p.stderr)
    else:
        colortext.warning(p.stdout)

def test_results(output_dir, PredictionSet):
    PredictionIDs = []
    results = get_results(PredictionSet)
    mutation_lists = {}
    for r in results:
        PredictionIDs.append(r['PredictionID'])
        mutation_lists[r['PredictionID']] = r['FlattenedMutations']
    RandomPredictionIDs = [PredictionIDs[random.randint(0, len(PredictionIDs) - 1)] for k in range(10)]
    RandomPredictionIDs = [54090L, 53875L, 54085L, 54079L, 54008L, 53853L, 53952L, 54056L, 53935L, 53893L]

    # Retrieve and unzip results
    if not(os.path.exists(output_dir)):
        os.mkdir(output_dir)
    for PredictionID in PredictionIDs:#RandomPredictionIDs:
        if not(os.path.exists(os.path.join(output_dir, str(PredictionID)))):
            colortext.message('Retrieving archive for Prediction %d.' % PredictionID)
            extract_data(ddG_connection, output_dir, PredictionID)

    # Get the sequences of the wildtype and mutant structures
    count = 0
    for PredictionID in PredictionIDs:#RandomPredictionIDs:
        wildtype_sequences = set()
        mutation_sequences = set()
        working_dir = os.path.join(os.path.join(output_dir, str(PredictionID)))
        for f in glob.glob(os.path.join(working_dir, '*.pdb')):
            if os.path.split(f)[1].startswith('mut_'):
                p = PDB.from_filepath(f)
                assert(len(p.atom_sequences) == 1)
                sequence = str(p.atom_sequences.values()[0])
                mutation_sequences.add(sequence)
            elif os.path.split(f)[1].startswith('repacked_wt_'):
                p = PDB.from_filepath(f)
                assert(len(p.atom_sequences) == 1)
                sequence = str(p.atom_sequences.values()[0])
                wildtype_sequences.add(sequence)

        assert(len(wildtype_sequences) == 1)
        assert(len(mutation_sequences) == 1)
        wildtype_sequence = wildtype_sequences.pop()
        mutation_sequence = mutation_sequences.pop()

        colortext.message('Prediction %d. Mutations: %s' % (PredictionID, mutation_lists[PredictionID]))
        assert(len(wildtype_sequence) == len(mutation_sequence))
        s = ''
        t = ''
        for x in range(len(wildtype_sequence)):
            if wildtype_sequence[x] != mutation_sequence[x]:
                s += colortext.make(wildtype_sequence[x], color="green")
                t += colortext.make(mutation_sequence[x], color="yellow")
            else:
                s += wildtype_sequence[x]
                t += mutation_sequence[x]
        print(s)
        print(t)

def create_pymol_session(output_filepath, download_dir, PredictionID, task_number, keep_files = True):
    '''Create a PyMOL session for a pair of structures.'''

    # Retrieve and unzip results
    if not(os.path.exists(download_dir)):
        os.mkdir(download_dir)
    working_dir = os.path.join(os.path.join(download_dir, str(PredictionID)))
    if not(os.path.exists(working_dir)) or not(os.path.exists(os.path.join(working_dir, 'repacked_wt_round_%d.pdb' % task_number))):
        extract_data(ddG_connection, download_dir, PredictionID)
    if not(os.path.exists(working_dir)) or not(os.path.exists(os.path.join(working_dir, 'repacked_wt_round_%d.pdb' % task_number))):
        raise Exception('Could not extract the models for task #%d of Prediction #%d.' % (task_number, PredictionID))

    # Retrieve the two structures corresponding to the task_number
    files = sorted(glob.glob(os.path.join(working_dir, '*_round_%d.pdb' % task_number)), reverse = True)
    assert(os.path.split(files[0])[1].startswith('repacked_wt_'))
    assert(os.path.split(files[1])[1].startswith('mut_'))

    # Creator the alignment object and write the PSE file
    chain_mapper = ScaffoldModelChainMapper.from_file_paths(files[0], files[1])

    # Remove the downloaded files
    if not keep_files:
        shutil.rmtree(download_dir)
    return chain_mapper.generate_pymol_session()

def write_pymol_session(output_filepath, download_dir, PredictionID, task_number, keep_files = True):
    PSE_file = create_pymol_session(output_filepath, download_dir, PredictionID, task_number, keep_files = keep_files)
    write_file(output_filepath, PSE_file[0], 'wb')

if False:

    # I have (hopefully) written these functions so that they can be run multiple times without consequences e.g.
    #  - if you accidentally add the same mutation string a second time, it will be handled gracefully;
    #  - if you run add_cluster_jobs multiple times, it won't force jobs to be re-run.

    # Here is an example usage:

    ### Computational stage ###

    # Step 1: Add a new PDB
    add_pdb_file('/kortemmelab/home/oconchus/BiosensorDesign/S9G10_best.pdb', 'S9G10_best')
    add_pdb_file('/kortemmelab/home/oconchus/BiosensorDesign/1SVX.pdb', '1SVX')

    # Step 2: Add a list of mutations for this PDB
    #
    # generate_all_combinations_of_mutations generates all possible combinations of the list of mutations which can be useful.

    mutation_set_1 = [
        ChainMutation('F',  '70', 'Y', Chain = 'A'),
        #ChainMutation('L',  '78', 'Y', Chain = 'A'),
        ChainMutation('W', '103', 'M', Chain = 'A'),
        ChainMutation('W', '111', 'K', Chain = 'A'),
        ChainMutation('E', '141', 'L', Chain = 'A'),
        ChainMutation('S', '142', 'A', Chain = 'A'),
        ChainMutation('L', '150', 'I', Chain = 'A'),
        ChainMutation('I', '153', 'V', Chain = 'A'),
    ]
    mutant_list = generate_all_combinations_of_mutations(mutation_set_1)
    colortext.message('Adding %d mutants.' % len(mutant_list))
    for mutant_mutations in mutant_list:
        add_mutant('S9G10_best', mutant_mutations)

    mutation_set_2 = [
        ChainMutation('W', '103', 'M', Chain = 'A'),
        ChainMutation('W', '111', 'K', Chain = 'A'),
        ChainMutation('W', '111', 'L', Chain = 'A'),
        ChainMutation('E', '141', 'L', Chain = 'A'),
        ChainMutation('L', '150', 'I', Chain = 'A'),
    ]
    for mutant_mutations in generate_all_combinations_of_mutations(mutation_set_2):
        add_mutant('S9G10_best', mutant_mutations)

    # Step 2: If you wanted to add a single mutant:
    add_mutant('S9G10_best', [ChainMutation('F',  '70', 'Y', Chain = 'A'), ChainMutation('W', '103', 'M', Chain = 'A')]) # pdb_id, chain, comma-separated list of mutations


    # Step 3:
    # This function will queue any mutations of 3SO8_BS1 that we added with add_mutations which have not previously run on
    # the cluster or which are not in the queue.
    # The second parameter is a 'prediction set'. This is just a label used to group a set of predictions together. It makes sense
    # to group all predictions using the same protocol together.
    # The third parameter is the name of the DDG protocol. All protocols currently in the database are variants of protocol 16
    # from the Kellogg, Leaver-Fay, Baker paper (doi:10.1002/prot.22921).

    add_cluster_jobs('S9G10_best', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)')
    #add_cluster_jobs('1SVX', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)')


    # Steps 2 and 3 can be repeated as often as you need however it is best to add as many mutations in step 2 first as this
    # will result in a better use of the cluster (larger array jobs).

    # todo: remove 1MJ0_BS1 used for testing

    ### Analysis stage ###

    # To see the results quickly, you can use the get_results function e.g.
    colortext.message('Retrieving results')
    get_results('FPP biosensor: protocol 16')

    # The analyze_results function will create a sorted graph of the results showing which sets of mutations are predicted to be more stable
    colortext.message('Analyzing results')
    analyze_results('FPP biosensor: protocol 16', 'L87Y_removed.png', 'kellogg', 'total')

    # This is a simple test function which prints out the sequences of the monomer wildtype sequence and mutant sequence,
    # highlighting where they differ in case there was a mistake in the pipeline.
    # First, the output files are downloaded and extracted to the directory specified in the first argument.
    colortext.message('Testing results')
    test_results('random_output_data', 'FPP biosensor: protocol 16')

    # This function creates a PyMOL session
    # The output files are downloaded and extracted to the directory specified in the second argument.
    # This function takes in a Prediction ID. These can be retrieved using the get_results function above.
    # The fourth argument is a task ID e.g. Protocol 16 generates 50 pairs of wildtype and mutant models by default, numbered
    # 1 to 50. This argument picks one of those pairs and creates a PyMOL session with aligned structures and with the mutations
    # highlighted.
    colortext.message('Creating PyMOL session')
    create_pymol_session('test.pse', 'random_output_data', 53858, 25)

#analyze_results('FPP biosensor: protocol 16', 'L87Y_removed_Noah.png', 'noah_8,0A', 'positional')
#analyze_results('FPP biosensor: protocol 16', 'L87Y_removed.png', 'kellogg', 'total')
#add_cluster_jobs('S9G10_best', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)')