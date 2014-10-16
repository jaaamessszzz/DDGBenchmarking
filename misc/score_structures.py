# A script to retrieve the structure scores from the prediction's standard output file and store them in the database

import sys
sys.path.insert(0, "../..")
sys.path.insert(0, "../ddglib")
import zipfile
import shutil
import os
import time
import pickle
import subprocess
import json
from tools import colortext
from tools.deprecated.rosettahelper import readBinaryFile, makeTemp755Directory, writeFile, readFileLines
from tools.bio.pdb import PDB
from tools.bio.basics import residue_type_3to1_map as aa1

import ddgdbapi
from dbapi import ddG as ddGInterface


def get_completed_prediction_sets(DDG_api):
    ddGdb = DDG_api.ddGDB
    completed_prediction_sets = []
    prediction_set_statuses = {}
    results = ddGdb.execute('SELECT PredictionSet, status, COUNT(ID) FROM Prediction WHERE PredictionSet LIKE "%ubiquitin%" GROUP BY PredictionSet, status')
    for r in results:
        prediction_set_statuses[r['PredictionSet']] = prediction_set_statuses.get(r['PredictionSet'], set())
        prediction_set_statuses[r['PredictionSet']].add(r['status'])
    for prediction_set, statuses in sorted(prediction_set_statuses.iteritems()):
        if len(statuses) == 1 and statuses.pop() == 'done':
            completed_prediction_sets.append(prediction_set)
    return completed_prediction_sets


def determine_structure_scores(DDG_api):
    ddGdb = DDG_api.ddGDB

    # Get the list of completed prediction set
    completed_prediction_sets = get_completed_prediction_sets(DDG_api)

    # For each completed prediction set, determine the structure scores
    for prediction_set in completed_prediction_sets:
        predictions = ddGdb.execute('SELECT ID, status FROM Prediction WHERE PredictionSet=%s AND StructureScores IS NULL', parameters=(prediction_set,))
        num_predictions = len(predictions)
        count = 1

        # Iterate over all completed Predictions with null StructureScores. For each Prediction, determine and store the structure scores
        for prediction in predictions:
            colortext.message('%s: %d of %d' % (prediction_set, count, num_predictions))
            assert(prediction['status'] == 'done')
            PredictionID = prediction['ID']

            # Get the ddg_monomer scores for each structure
            grouped_scores = json.dumps(DDG_api.get_ddg_monomer_scores_per_structure(PredictionID))

            # Store this mapping as a JSON string in the database
            ddGdb.execute('UPDATE Prediction SET StructureScores=%s WHERE ID=%s', parameters=(grouped_scores, PredictionID))

            count += 1
            break
        break


if __name__ == '__main__':
    DDG_api = ddGInterface()
    #ddGdb = ddgdbapi.ddGDatabase()
    #ddGPredictiondb = ddgdbapi.ddGPredictionDataDatabase()
    determine_structure_scores(DDG_api)
