import sys
import os
import pprint
#sys.path.insert(0, "..")
sys.path.insert(0, "../ddglib")
sys.path.insert(0, "../..")
sys.path.insert(0, "..")

from ddglib import analysis
from tools import colortext
from tools.fs.fsio import read_file, get_file_lines, write_file
import json

filetype = "pdf"
prediction_set_for_analysis = 'Protocol_16_r57471'
prediction_set_for_analysis = 'RosCon2013_P16_score12prime'

colortext.message("*** %s ***" % prediction_set_for_analysis)

from ddglib import dbapi, ddgdbapi
from tools.bio import pdb



def dump_data():
    ddG_connection = dbapi.ddG()
    ddGdb = ddgdbapi.ddGDatabase()

    prediction_set = 'Protocol_16_r57471'
    userdata_set = 'AllValidPGPK'

    cached_pdb_details = json.loads(read_file('cached_pdb_details.json'))
    analysis_breakdown = ddG_connection.get_predictionset_data(prediction_set, userdata_set, cached_pdb_details = cached_pdb_details)

    test_data = dict(
        amino_acids = analysis_breakdown.amino_acids,
        pdb_details = analysis_breakdown.pdb_details,
        predictions = analysis_breakdown.predictions,
        analysis_datasets = analysis_breakdown.analysis_datasets,
    )
    write_file('example_analysis_input.json', json.dumps(test_data))


def load_data():
    d = json.loads(read_file('example_analysis_input.json'))
    return dbapi.AnalysisBreakdown(d['amino_acids'], d['pdb_details'], d['predictions'], d['analysis_datasets'])


if __name__ == '__main__':
    import json

    # Retrieve the data
    analysis_breakdown = load_data()

    for analysis_subset in analysis_breakdown.analysis_datasets.keys():
        #if analysis_subset != 'Kellogg':
        #    continue
        if analysis_subset == 'ProTherm':
            continue
        for scoring_method in ['Kellogg', 'Noah']:
            #analysis_breakdown.analyze_subset_main(analysis_subset, scoring_method)
            #analysis_breakdown.analyze_subset_GP(analysis_subset, scoring_method)
            #analysis_breakdown.analyze_subset_by_mutation_size(analysis_subset, scoring_method)
            #analysis_breakdown.analyze_subset_by_secondary_structure(analysis_subset, scoring_method)
            analysis_breakdown.analyze_subset_by_binned_chain_length(analysis_subset, scoring_method, num_bins = 9)

            #analysis_breakdown.analyze_subset_by_aromaticity(analysis_subset, scoring_method)
            #analysis_breakdown.analyze_subset_by_exposure(analysis_subset, scoring_method, cut_off=0.4)
            #analysis_breakdown.analyze_subset_by_specific_resolutions(analysis_subset, scoring_method)
            #analysis_breakdown.analyze_subset_by_binned_resolutions(analysis_subset, scoring_method)
        #break
