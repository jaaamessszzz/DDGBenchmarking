import sys
import datetime
sys.path.insert(0, "..")

from tools import colortext
from tools.bio.pdb import PDB
from tools.bio.basics import Mutation
from tools.fs.fsio import read_file
from ddglib import ddgdbapi, dbapi

DDGdb = ddgdbapi.ddGDatabase()
ddG_connection = dbapi.ddG()

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
    existing_pdb = DDGdb.execute('SELECT ID FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))

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
    ExperimentIDs = set([r['ID'] for r in DDGdb.execute('SELECT ID FROM Experiment WHERE PDBFileID=%s', parameters=(pdb_ID,))])
    ExperimentIDsInPredictionSet = set([r['ExperimentID'] for r in DDGdb.execute('SELECT ExperimentID FROM Prediction WHERE PredictionSet=%s', parameters=(PredictionSet,))])
    experiment_IDs_to_add = sorted(ExperimentIDs.difference(ExperimentIDsInPredictionSet))

    if experiment_IDs_to_add:
        colortext.printf("\nAdding %d jobs to the cluster queue." % len(experiment_IDs_to_add), "lightgreen")

        for experiment_ID in experiment_IDs_to_add:
            colortext.write('.', "lightgreen")
            ddG_connection.addPrediction(experiment_ID, None, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True)
    else:
        colortext.printf("\nAll jobs are already in the queue or have been run.", "lightgreen")
    print('')


# I have (hopefully) written these functions so that they can be run multiple times without consequences e.g.
#  - if you accidentally add the same mutation string a second time, it will be handled gracefully;
#  - if you run add_cluster_jobs multiple times, it won't force jobs to be re-run.

# Here is an example usage:

# Step 1: Add a new PDB
add_pdb_file('/kortemmelab/home/oconchus/BiosensorDesign/S9G10_best.pdb', 'S9G10_best')

# Step 2: Add a list of mutations for this PDB
add_mutations_from_string('S9G10_best', 'A', all_mutations)
add_mutations('S9G10_best', 'A', 'F70Y, W111K,') # pdb_id, chain, comma-separated list of mutations

# Step 3:
# This function will queue any mutations of 3SO8_BS1 that we added with add_mutations which have not previously run on
# the cluster or which are not in the queue.
# The second parameter is a 'prediction set'. This is just a label used to group a set of predictions together. It makes sense
# to group all predictions using the same protocol together.
# The third parameter is the name of the DDG protocol. All protocols currently in the database are variants of protocol 16
# from the Kellogg, Leaver-Fay, Baker paper (doi:10.1002/prot.22921).

add_cluster_jobs('S9G10_best', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)')

# Steps 2 and 3 can be repeated as often as you need however it is best to add as many mutations in step 2 first as this
# will result in a better use of the cluster (larger array jobs).

# todo: remove 1MJ0_BS1 used for testing