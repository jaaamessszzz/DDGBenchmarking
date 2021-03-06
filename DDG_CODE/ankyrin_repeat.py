import sys
sys.path.insert(0, "..")


from klab import colortext
from klab.bio.pdb import PDB
from klab.bio.basics import Mutation, ChainMutation, generate_all_combinations_of_mutations
from klab.bio.alignment import ScaffoldModelChainMapper
from klab.fs.fsio import read_file, write_file

from ddglib import ddgdbapi, db_api

if __name__ == '__main__':
    DDGdb = ddgdbapi.ddGDatabase()
    ddG_connection = db_api.ddG()

all_wildtype_mutations = '''
# All mutations are on chain A
A83P, I86L, S87A
L95I, I98V
A83P, I86L, S87A, L95I, I98V
'''

if False:

    # I have (hopefully) written these functions so that they can be run multiple times without consequences e.g.
    #  - if you accidentally add the same mutation string a second time, it will be handled gracefully;
    #  - if you run add_predictions_by_pdb_id multiple times, it won't force jobs to be re-run.

    # Here is an example usage:

    ### Computational stage ###

    # Step 1: Add a new PDB
    ddG_connection.add_pdb_file('/kortemmelab/home/oconchus/BiosensorDesign/S9G10_best.pdb', 'S9G10_best')
    ddG_connection.add_pdb_file('/kortemmelab/home/oconchus/BiosensorDesign/1SVX.pdb', '1SVX')

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
        ddG_connection.add_mutant('S9G10_best', mutant_mutations)

    mutation_set_2 = [
        ChainMutation('W', '103', 'M', Chain = 'A'),
        ChainMutation('W', '111', 'K', Chain = 'A'),
        ChainMutation('W', '111', 'L', Chain = 'A'),
        ChainMutation('E', '141', 'L', Chain = 'A'),
        ChainMutation('L', '150', 'I', Chain = 'A'),
    ]
    for mutant_mutations in generate_all_combinations_of_mutations(mutation_set_2):
        ddG_connection.add_mutant('S9G10_best', mutant_mutations)

    # Step 2: If you wanted to add a single mutant:
    ddG_connection.add_mutant('S9G10_best', [ChainMutation('F',  '70', 'Y', Chain = 'A'), ChainMutation('W', '103', 'M', Chain = 'A')]) # pdb_id, chain, comma-separated list of mutations


    # Step 3:
    # This function will queue any mutations of 3SO8_BS1 that we added with add_mutations which have not previously run on
    # the cluster or which are not in the queue.
    # The second parameter is a 'prediction set'. This is just a label used to group a set of predictions together. It makes sense
    # to group all predictions using the same protocol together.
    # The third parameter is the name of the DDG protocol. All protocols currently in the database are variants of protocol 16
    # from the Kellogg, Leaver-Fay, Baker paper (doi:10.1002/prot.22921).

    ddG_connection.add_predictions_by_pdb_id('S9G10_best', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)', priority = 9)
    #ddG_connection.add_predictions_by_pdb_id('1SVX', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)', priority = 9)


    # Steps 2 and 3 can be repeated as often as you need however it is best to add as many mutations in step 2 first as this
    # will result in a better use of the cluster (larger array jobs).

    # todo: remove 1MJ0_BS1 used for testing

    ### Analysis stage ###

    # To see the results quickly, you can use the get_flattened_prediction_results function e.g.
    colortext.message('Retrieving results')
    ddG_connection.get_flattened_prediction_results('FPP biosensor: protocol 16')

    # The create_abacus_graph_for_a_single_structure function will create a sorted graph of the results showing which sets of mutations are predicted to be more stable
    colortext.message('Analyzing results')
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'L87Y_removed.png')

    # This is a simple test function which prints out the sequences of the monomer wildtype sequence and mutant sequence,
    # highlighting where they differ in case there was a mistake in the pipeline.
    # First, the output files are downloaded and extracted to the directory specified in the first argument.
    colortext.message('Testing results')
    ddG_connection.test_results('random_output_data', 'FPP biosensor: protocol 16')

    # This function creates a PyMOL session
    # The output files are downloaded and extracted to the directory specified in the second argument.
    # This function takes in a Prediction ID. These can be retrieved using the get_flattened_prediction_results function above.
    # The fourth argument is a task ID e.g. Protocol 16 generates 50 pairs of wildtype and mutant models by default, numbered
    # 1 to 50. This argument picks one of those pairs and creates a PyMOL session with aligned structures and with the mutations
    # highlighted.
    colortext.message('Creating PyMOL session')
    ddG_connection.create_pymol_session('test.pse', 'random_output_data', 53858, 25)

#create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'noah_8,0A', 'positional', graph_filename = 'L87Y_removed_Noah.png')
#create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'L87Y_removed.png')
#ddG_connection.add_predictions_by_pdb_id('S9G10_best', 'FPP biosensor: protocol 16', 'Protocol16 3.5.1 (talaris2013sc)', priority = 9)
#ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'L87Y_removed_new.png')

def test_abacus_graph():
    '''This function can be deleted. It was added to test the abacus graph with different numbers of datapoints.'''
    import os
    import json
    if not(os.path.exists('results_cache.txt')):
        results = ddG_connection.get_flattened_prediction_results('FPP biosensor: protocol 16')
        for r in results:
            r['TimeTaken'] = r['TimeTaken'].total_seconds() # timedelta objects are not JSON serializable
        write_file('results_cache.txt', json.dumps(results), 'w')
    results = json.loads(read_file('results_cache.txt'))

    try:
        ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_3.png', cached_results = results, num_datapoints = 3)
    except:
        pass
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_5.png', cached_results = results, num_datapoints = 5)
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_8.png', cached_results = results, num_datapoints = 8)
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_10.png', cached_results = results, num_datapoints = 10)
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_12.png', cached_results = results, num_datapoints = 12)
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_20.png', cached_results = results, num_datapoints = 20)
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_50.png', cached_results = results, num_datapoints = 50)
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_127.png', cached_results = results, num_datapoints = 127)
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'test_255.png', cached_results = results, num_datapoints = 255)

#test_abacus_graph()

#p = PDB.from_filepath('/kortemmelab/home/oconchus/BiosensorDesign/S9G10_best.pdb')

#print(p.lines)
#b_seq = p.atom_sequences['B']
#print(' '.join([res.get_residue_id().replace(' ', '') for id, res in b_seq]))

#ddG_connection.add_predictions_by_pdb_id('S9G10_best', 'FPP biosensor: p16 (complex) full', 'Protocol16 3.5.1 (FPP complex) (talaris2013sc)', priority = 9, KeepHETATMLines = True, strip_other_chains = False, status = 'halted')
