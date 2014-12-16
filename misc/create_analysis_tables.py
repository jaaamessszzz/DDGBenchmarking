import sys
import os
#sys.path.insert(0, "..")
sys.path.insert(0, "../ddglib")
sys.path.insert(0, "../..")
sys.path.insert(0, "..")

from ddglib import analysis
from tools import colortext
from tools.fs.fsio import read_file, get_file_lines

filetype = "pdf"
prediction_set_for_analysis = 'Protocol_16_r57471'
prediction_set_for_analysis = 'RosCon2013_P16_score12prime'

colortext.message("*** %s ***" % prediction_set_for_analysis)

from ddglib import dbapi, ddgdbapi
from tools.bio import pdb


if __name__ == '__main__':
    import json
    cached_pdb_details = json.loads(read_file('cached_pdb_details.json'))
    ddG_connection = dbapi.ddG()
    ddGdb = ddgdbapi.ddGDatabase()

    # Retrieve the data
    subset = 'Guerois'
    subset = 'AlaScan-GPK'
    data = ddG_connection.get_predictionset_data('Protocol_16_r57471', cached_pdb_details = cached_pdb_details)
    predictions = data['predictions']
    analysis_dataset = data['analysis_datasets'][subset]

    # Analyze the data
    from tools.stats.misc import get_xy_dataset_correlations
    xvalues = []
    yvalues = []
    for prediction_id, details in sorted(predictions.iteritems()):
        if prediction_id in analysis_dataset:
            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
            predicted_score = details['Kellogg']
            if predicted_score != None:
                xvalues.append(ExperimentalDDG)
                yvalues.append(predicted_score)

    colortext.message('Analyzing %d values for dataset %s.' % (len(xvalues), subset))
    print(get_xy_dataset_correlations(xvalues, yvalues))


