import sys
import os
import datetime
import string
import random
import pickle
import glob
import shutil
sys.path.insert(0, "..")

import numpy as np
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt

from tools import colortext
from tools.bio.pdb import PDB
from tools.bio.basics import Mutation
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

all_mutations = '''
# All mutations are on chain A
F70Y,
W103M,
L78Y,
W111K,
E141L,
S142A,
L150I,
I153V,
F70Y, W103M,
F70Y, L78Y,
F70Y, W111K,
F70Y, E141L,
F70Y, S142A,
F70Y, L150I,
F70Y, I153V,
W103M, L78Y,
W103M, W111K,
W103M, E141L,
W103M, S142A,
W103M, L150I,
W103M, I153V,
L78Y, W111K,
L78Y, E141L,
L78Y, S142A,
L78Y, L150I,
L78Y, I153V,
W111K, E141L,
W111K, S142A,
W111K, L150I,
W111K, I153V,
E141L, S142A,
E141L, L150I,
E141L, I153V,
S142A, L150I,
S142A, I153V,
L150I, I153V,
F70Y, W103M, L78Y,
F70Y, W103M, W111K,
F70Y, W103M, E141L,
F70Y, W103M, S142A,
F70Y, W103M, L150I,
F70Y, W103M, I153V,
F70Y, L78Y, W111K,
F70Y, L78Y, E141L,
F70Y, L78Y, S142A,
F70Y, L78Y, L150I,
F70Y, L78Y, I153V,
F70Y, W111K, E141L,
F70Y, W111K, S142A,
F70Y, W111K, L150I,
F70Y, W111K, I153V,
F70Y, E141L, S142A,
F70Y, E141L, L150I,
F70Y, E141L, I153V,
F70Y, S142A, L150I,
F70Y, S142A, I153V,
F70Y, L150I, I153V,
W103M, L78Y, W111K,
W103M, L78Y, E141L,
W103M, L78Y, S142A,
W103M, L78Y, L150I,
W103M, L78Y, I153V,
W103M, W111K, E141L,
W103M, W111K, S142A,
W103M, W111K, L150I,
W103M, W111K, I153V,
W103M, E141L, S142A,
W103M, E141L, L150I,
W103M, E141L, I153V,
W103M, S142A, L150I,
W103M, S142A, I153V,
W103M, L150I, I153V,
L78Y, W111K, E141L,
L78Y, W111K, S142A,
L78Y, W111K, L150I,
L78Y, W111K, I153V,
L78Y, E141L, S142A,
L78Y, E141L, L150I,
L78Y, E141L, I153V,
L78Y, S142A, L150I,
L78Y, S142A, I153V,
L78Y, L150I, I153V,
W111K, E141L, S142A,
W111K, E141L, L150I,
W111K, E141L, I153V,
W111K, S142A, L150I,
W111K, S142A, I153V,
W111K, L150I, I153V,
E141L, S142A, L150I,
E141L, S142A, I153V,
E141L, L150I, I153V,
S142A, L150I, I153V,
F70Y, W103M, L78Y, W111K,
F70Y, W103M, L78Y, E141L,
F70Y, W103M, L78Y, S142A,
F70Y, W103M, L78Y, L150I,
F70Y, W103M, L78Y, I153V,
F70Y, W103M, W111K, E141L,
F70Y, W103M, W111K, S142A,
F70Y, W103M, W111K, L150I,
F70Y, W103M, W111K, I153V,
F70Y, W103M, E141L, S142A,
F70Y, W103M, E141L, L150I,
F70Y, W103M, E141L, I153V,
F70Y, W103M, S142A, L150I,
F70Y, W103M, S142A, I153V,
F70Y, W103M, L150I, I153V,
F70Y, L78Y, W111K, E141L,
F70Y, L78Y, W111K, S142A,
F70Y, L78Y, W111K, L150I,
F70Y, L78Y, W111K, I153V,
F70Y, L78Y, E141L, S142A,
F70Y, L78Y, E141L, L150I,
F70Y, L78Y, E141L, I153V,
F70Y, L78Y, S142A, L150I,
F70Y, L78Y, S142A, I153V,
F70Y, L78Y, L150I, I153V,
F70Y, W111K, E141L, S142A,
F70Y, W111K, E141L, L150I,
F70Y, W111K, E141L, I153V,
F70Y, W111K, S142A, L150I,
F70Y, W111K, S142A, I153V,
F70Y, W111K, L150I, I153V,
F70Y, E141L, S142A, L150I,
F70Y, E141L, S142A, I153V,
F70Y, E141L, L150I, I153V,
F70Y, S142A, L150I, I153V,
W103M, L78Y, W111K, E141L,
W103M, L78Y, W111K, S142A,
W103M, L78Y, W111K, L150I,
W103M, L78Y, W111K, I153V,
W103M, L78Y, E141L, S142A,
W103M, L78Y, E141L, L150I,
W103M, L78Y, E141L, I153V,
W103M, L78Y, S142A, L150I,
W103M, L78Y, S142A, I153V,
W103M, L78Y, L150I, I153V,
W103M, W111K, E141L, S142A,
W103M, W111K, E141L, L150I,
W103M, W111K, E141L, I153V,
W103M, W111K, S142A, L150I,
W103M, W111K, S142A, I153V,
W103M, W111K, L150I, I153V,
W103M, E141L, S142A, L150I,
W103M, E141L, S142A, I153V,
W103M, E141L, L150I, I153V,
W103M, S142A, L150I, I153V,
L78Y, W111K, E141L, S142A,
L78Y, W111K, E141L, L150I,
L78Y, W111K, E141L, I153V,
L78Y, W111K, S142A, L150I,
L78Y, W111K, S142A, I153V,
L78Y, W111K, L150I, I153V,
L78Y, E141L, S142A, L150I,
L78Y, E141L, S142A, I153V,
L78Y, E141L, L150I, I153V,
L78Y, S142A, L150I, I153V,
W111K, E141L, S142A, L150I,
W111K, E141L, S142A, I153V,
W111K, E141L, L150I, I153V,
W111K, S142A, L150I, I153V,
E141L, S142A, L150I, I153V,
F70Y, W103M, L78Y, W111K, E141L,
F70Y, W103M, L78Y, W111K, S142A,
F70Y, W103M, L78Y, W111K, L150I,
F70Y, W103M, L78Y, W111K, I153V,
F70Y, W103M, L78Y, E141L, S142A,
F70Y, W103M, L78Y, E141L, L150I,
F70Y, W103M, L78Y, E141L, I153V,
F70Y, W103M, L78Y, S142A, L150I,
F70Y, W103M, L78Y, S142A, I153V,
F70Y, W103M, L78Y, L150I, I153V,
F70Y, W103M, W111K, E141L, S142A,
F70Y, W103M, W111K, E141L, L150I,
F70Y, W103M, W111K, E141L, I153V,
F70Y, W103M, W111K, S142A, L150I,
F70Y, W103M, W111K, S142A, I153V,
F70Y, W103M, W111K, L150I, I153V,
F70Y, W103M, E141L, S142A, L150I,
F70Y, W103M, E141L, S142A, I153V,
F70Y, W103M, E141L, L150I, I153V,
F70Y, W103M, S142A, L150I, I153V,
F70Y, L78Y, W111K, E141L, S142A,
F70Y, L78Y, W111K, E141L, L150I,
F70Y, L78Y, W111K, E141L, I153V,
F70Y, L78Y, W111K, S142A, L150I,
F70Y, L78Y, W111K, S142A, I153V,
F70Y, L78Y, W111K, L150I, I153V,
F70Y, L78Y, E141L, S142A, L150I,
F70Y, L78Y, E141L, S142A, I153V,
F70Y, L78Y, E141L, L150I, I153V,
F70Y, L78Y, S142A, L150I, I153V,
F70Y, W111K, E141L, S142A, L150I,
F70Y, W111K, E141L, S142A, I153V,
F70Y, W111K, E141L, L150I, I153V,
F70Y, W111K, S142A, L150I, I153V,
F70Y, E141L, S142A, L150I, I153V,
W103M, L78Y, W111K, E141L, S142A,
W103M, L78Y, W111K, E141L, L150I,
W103M, L78Y, W111K, E141L, I153V,
W103M, L78Y, W111K, S142A, L150I,
W103M, L78Y, W111K, S142A, I153V,
W103M, L78Y, W111K, L150I, I153V,
W103M, L78Y, E141L, S142A, L150I,
W103M, L78Y, E141L, S142A, I153V,
W103M, L78Y, E141L, L150I, I153V,
W103M, L78Y, S142A, L150I, I153V,
W103M, W111K, E141L, S142A, L150I,
W103M, W111K, E141L, S142A, I153V,
W103M, W111K, E141L, L150I, I153V,
W103M, W111K, S142A, L150I, I153V,
W103M, E141L, S142A, L150I, I153V,
L78Y, W111K, E141L, S142A, L150I,
L78Y, W111K, E141L, S142A, I153V,
L78Y, W111K, E141L, L150I, I153V,
L78Y, W111K, S142A, L150I, I153V,
L78Y, E141L, S142A, L150I, I153V,
W111K, E141L, S142A, L150I, I153V,
F70Y, W103M, L78Y, W111K, E141L, S142A,
F70Y, W103M, L78Y, W111K, E141L, L150I,
F70Y, W103M, L78Y, W111K, E141L, I153V,
F70Y, W103M, L78Y, W111K, S142A, L150I,
F70Y, W103M, L78Y, W111K, S142A, I153V,
F70Y, W103M, L78Y, W111K, L150I, I153V,
F70Y, W103M, L78Y, E141L, S142A, L150I,
F70Y, W103M, L78Y, E141L, S142A, I153V,
F70Y, W103M, L78Y, E141L, L150I, I153V,
F70Y, W103M, L78Y, S142A, L150I, I153V,
F70Y, W103M, W111K, E141L, S142A, L150I,
F70Y, W103M, W111K, E141L, S142A, I153V,
F70Y, W103M, W111K, E141L, L150I, I153V,
F70Y, W103M, W111K, S142A, L150I, I153V,
F70Y, W103M, E141L, S142A, L150I, I153V,
F70Y, L78Y, W111K, E141L, S142A, L150I,
F70Y, L78Y, W111K, E141L, S142A, I153V,
F70Y, L78Y, W111K, E141L, L150I, I153V,
F70Y, L78Y, W111K, S142A, L150I, I153V,
F70Y, L78Y, E141L, S142A, L150I, I153V,
F70Y, W111K, E141L, S142A, L150I, I153V,
W103M, L78Y, W111K, E141L, S142A, L150I,
W103M, L78Y, W111K, E141L, S142A, I153V,
W103M, L78Y, W111K, E141L, L150I, I153V,
W103M, L78Y, W111K, S142A, L150I, I153V,
W103M, L78Y, E141L, S142A, L150I, I153V,
W103M, W111K, E141L, S142A, L150I, I153V,
L78Y, W111K, E141L, S142A, L150I, I153V,
F70Y, W103M, L78Y, W111K, E141L, S142A, L150I,
F70Y, W103M, L78Y, W111K, E141L, S142A, I153V,
F70Y, W103M, L78Y, W111K, E141L, L150I, I153V,
F70Y, W103M, L78Y, W111K, S142A, L150I, I153V,
F70Y, W103M, L78Y, E141L, S142A, L150I, I153V,
F70Y, W103M, W111K, E141L, S142A, L150I, I153V,
F70Y, L78Y, W111K, E141L, S142A, L150I, I153V,
W103M, L78Y, W111K, E141L, S142A, L150I, I153V,
F70Y, W103M, L78Y, W111K, E141L, S142A, L150I, I153V,
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


def add_mutations_from_file(pdb_ID, chain, filepath):
    '''filepath should have the following format:
         - each non-empty line which does not begin with '#' should be a comma-separated list of mutations in the format
           wtAA_resID_mutantAA where resID can contain a insertion code e.g. P50A. Each line corresponds to one 'experiment'
           so a line with one mutation corresponds to a single mutant, a line with two mutations corresponds to a double mutant etc.
         - empty lines are lines beginning with '#' are ignored.'''
    add_mutations_from_string(pdb_ID, chain, read_file(filepath))


def add_mutations_from_string(pdb_ID, chain, str):
    '''str should have the following format:
         - each non-empty line which does not begin with '#' should be a comma-separated list of mutations in the format
           wtAA_resID_mutantAA where resID can contain a insertion code e.g. P50A. Each line corresponds to one 'experiment'
           so a line with one mutation corresponds to a single mutant, a line with two mutations corresponds to a double mutant etc.
         - empty lines are lines beginning with '#' are ignored.'''

    lines = []
    for l in str.split('\n'):
        l = l.strip()
        if l and not(l.startswith('#')):
            lines.append(l)
    for l in lines:
        add_mutations(pdb_ID, chain, l)


def add_mutations(pdb_ID, chain, mutation_string):
    '''Use this function to add one set of mutations (i.e. corresponding to one mutant) to the database.

       mutation_string here is a comma-separated list of mutations in the format wtAA_resID_mutantAA where resID can contain a insertion code
       (even though the model is from Rosetta in this case so there will be no insertion code).
    '''
    mutations = []
    for mutation_string in [m.strip() for m in mutation_string.split(',') if m.strip()]:
        mutations.append(Mutation(mutation_string[0], mutation_string[1:-1], mutation_string[-1], chain))

    colortext.warning("Adding mutation: %s." % ', '.join(map(str, mutations)))
    ddG_connection.createDummyExperiment_ankyrin_repeat(pdb_ID, mutations, chain)


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
SELECT Prediction.ID AS PredictionID, Prediction.ExperimentID, Experiment.PDBFileID, ExperimentMutations.FlattenedMutations, Prediction.ddG, TIMEDIFF(Prediction.EndDate, Prediction.StartDate) AS TimeTaken FROM Prediction INNER JOIN
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

def analyze_results(PredictionSet, graph_filename):
    results = get_results(PredictionSet)
    sorted_results = {}
    for r in results:
        sorted_results[(pickle.loads(r['ddG'])['data']['kellogg']['total']['ddG'], r['ExperimentID'])] = r
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
        #if r['FlattenedMutations'].find('A L78Y') == -1:
            #print('%f, %s' % (k[0], r['FlattenedMutations']))
            #count += 1

    data = []
    pruned_data = []
    for k, r in sorted(sorted_results.iteritems()):
        line = []
        for m in sorted(set_of_mutations):
            if r['FlattenedMutations'].find(m[1]) != -1:
                line.append(1)
            else:
                line.append(0)
        data.append((pickle.loads(r['ddG'])['data']['kellogg']['total']['ddG'], line))
        if r['FlattenedMutations'].find('A L78Y') == -1:
            pruned_data.append((pickle.loads(r['ddG'])['data']['kellogg']['total']['ddG'], line))

    labels = [m[1] for m in sorted(set_of_mutations)]

    #hinton('all_mutants.png', labels, data)
    create_graph(graph_filename, labels, pruned_data)

def create_graph(filename, labels, data, max_weight=None, ax=None):
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
        plt.text(x, -30, l.split()[1], fontdict=None, withdash=True, fontsize=9)
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
    plt.title(r'$\Delta\Delta$G predictions for FPP biosensor mutants', fontdict=None)

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
    #print(os.path.join(working_dir, '*.pdb'))
    # Retrieve and unzip results

    if not(os.path.exists(download_dir)):
        os.mkdir(download_dir)

    working_dir = os.path.join(os.path.join(download_dir, str(PredictionID)))
    if not(os.path.exists(working_dir)) or not(os.path.exists(os.path.join(working_dir, 'repacked_wt_round_%d.pdb' % task_number))):
        extract_data(ddG_connection, download_dir, PredictionID)
    if not(os.path.exists(working_dir)) or not(os.path.exists(os.path.join(working_dir, 'repacked_wt_round_%d.pdb' % task_number))):
        raise Exception('Could not extract the models for task #%d of Prediction #%d.' % (task_number, PredictionID))

    files = sorted(glob.glob(os.path.join(working_dir, '*_round_%d.pdb' % task_number)), reverse = True)
    assert(os.path.split(files[0])[1].startswith('repacked_wt_'))
    assert(os.path.split(files[1])[1].startswith('mut_'))

    chain_mapper = ScaffoldModelChainMapper.from_file_paths(files[0], files[1])
    PSE_file = chain_mapper.generate_pymol_session()

    write_file(output_filepath, PSE_file[0], 'wb')

    if not keep_files:
        shutil.rmtree(download_dir)

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
    add_mutations_from_string('S9G10_best', 'A', all_mutations)
    add_mutations_from_string('1SVX', 'A', all_wildtype_mutations)
    # Step 2: If you wanted to add a single mutation:
    add_mutations('S9G10_best', 'A', 'F70Y, W111K,') # pdb_id, chain, comma-separated list of mutations

    # Step 3:
    # This function will queue any mutations of 3SO8_BS1 that we added with add_mutations which have not previously run on
    # the cluster or which are not in the queue.
    # The second parameter is a 'prediction set'. This is just a label used to group a set of predictions together. It makes sense
    # to group all predictions using the same protocol together.
    # The third parameter is the name of the DDG protocol. All protocols currently in the database are variants of protocol 16
    # from the Kellogg, Leaver-Fay, Baker paper (doi:10.1002/prot.22921).

    #add_cluster_jobs('S9G10_best', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)')
    add_cluster_jobs('1SVX', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)')

    # Steps 2 and 3 can be repeated as often as you need however it is best to add as many mutations in step 2 first as this
    # will result in a better use of the cluster (larger array jobs).

    # todo: remove 1MJ0_BS1 used for testing

    ### Analysis stage ###

    # To see the results quickly, you can use the get_results function e.g.
    get_results('FPP biosensor: protocol 16')

    # The analyze_results function will create a sorted graph of the results showing which sets of mutations are predicted to be more stable
    analyze_results('FPP biosensor: protocol 16', 'L87Y_removed.png')

    # This is a simple test function which prints out the sequences of the monomer wildtype sequence and mutant sequence,
    # highlighting where they differ in case there was a mistake in the pipeline.
    # First, the output files are downloaded and extracted to the directory specified in the first argument.
    test_results('random_output_data', 'FPP biosensor: protocol 16')

    # This function creates a PyMOL session
    # The output files are downloaded and extracted to the directory specified in the second argument.
    # This function takes in a Prediction ID. These can be retrieved using the get_results function above.
    # The fourth argument is a task ID e.g. Protocol 16 generates 50 pairs of wildtype and mutant models by default, numbered
    # 1 to 50. This argument picks one of those pairs and creates a PyMOL session with aligned structures and with the mutations
    # highlighted.
    create_pymol_session('test.pse', 'random_output_data', 53858, 25)

test_results('random_output_data', 'FPP biosensor: protocol 16')
