#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

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
import pprint
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
    return ['RosCon2013_P16_talaris2013sc']
    return completed_prediction_sets


def determine_structure_scores(DDG_api, skip_if_we_have_pairs = 50):
    pp = pprint.PrettyPrinter(indent=4)

    ddGdb = DDG_api.ddGDB
    ddGdb_utf = ddgdbapi.ddGDatabase(use_utf = True)
    # Get the list of completed prediction set
    completed_prediction_sets = get_completed_prediction_sets(DDG_api)
    print(completed_prediction_sets)

    # Create the mapping from the old score types to the ScoreMethod record IDs
    ScoreMethodMap = {}
    results = ddGdb_utf.execute('SELECT * FROM ScoreMethod')
    for r in results:
        if r['MethodName'] == 'Global' and r['MethodType'] == 'Protocol 16':
            ScoreMethodMap[("kellogg", "total")] = r['ID']
        if r['Authors'] == 'Noah Ollikainen':
            if r['MethodName'] == 'Local' and r['MethodType'] == 'Position' and r['Parameters'] == u'8Å radius':
                ScoreMethodMap[("noah_8,0A", "positional")] = r['ID']
            if r['MethodName'] == 'Local' and r['MethodType'] == 'Position (2-body)' and r['Parameters'] == u'8Å radius':
                ScoreMethodMap[("noah_8,0A", "positional_twoscore")] = r['ID']
            if r['MethodName'] == 'Global' and r['MethodType'] == 'By residue' and r['Parameters'] == u'8Å radius':
                ScoreMethodMap[("noah_8,0A", "total")] = r['ID']

    # For each completed prediction set, determine the structure scores
    for prediction_set in completed_prediction_sets:
        #if prediction_set not in ['Ubiquitin scan: UQ_con_yeast p16']:
        #    continue

        predictions = ddGdb.execute('SELECT ID, ddG, Scores, status, ScoreVersion FROM Prediction WHERE PredictionSet=%s ORDER BY ID', parameters=(prediction_set,))
        num_predictions = len(predictions)

        # Pass #1: Iterate over all Predictions and make sure that they gave completed and contain all the scores we expect
        colortext.message('Prediction set: %s' % prediction_set)
        colortext.warning('Checking that all data exists...')
        for prediction in predictions:
            #assert(prediction['status'] == 'done')
            PredictionID = prediction['ID']
            if PredictionID != 72856:
                continue
            global_scores = pickle.loads(prediction['ddG'])
            assert(global_scores)
            assert(prediction['ScoreVersion'] == 0.23)
            if not prediction['Scores']:
                raise Exception("This prediction needs to be scored with Noah's method.")

            gs2 = json.loads(prediction['Scores'])
            if True not in set([k.find('noah') != -1 for k in gs2['data'].keys()]):
                raise Exception("This prediction needs to be scored with Noah's method.")
            assert (gs2['data']['kellogg'] == global_scores['data']['kellogg'])

        # Pass #2: Iterate over all completed Predictions with null StructureScores.
        # For each Prediction, determine and store the structure scores
        count = 0
        for prediction in predictions:

            count += 1
            PredictionID = prediction['ID']
            colortext.message('%s: %d of %d (Prediction #%d)' % (prediction_set, count, num_predictions, PredictionID))

            #if PredictionID != 72856:
            #if PredictionID < 73045: continue
            if prediction['status'] == 'failed':
                colortext.error('Skipping failed prediction %d.' % PredictionID)
                continue
            if prediction['status'] == 'queued':
                colortext.warning('Skipping queued prediction %d.' % PredictionID)
                continue
            if prediction['status'] == 'postponed':
                colortext.printf('Skipping postponed prediction %d.' % PredictionID, 'cyan')
                continue

            # Store the ensemble scores
            try:
                global_scores = json.loads(prediction['Scores'])['data']
            except:
                raise colortext.Exception("Failed reading the Scores field's JSON object. The Prediction Status is %(status)s. The Scores field is: '%(Scores)s'." % prediction)
            for score_type, inner_data in global_scores.iteritems():
                for inner_score_type, data in inner_data.iteritems():
                    components = {}
                    if score_type == 'kellogg' and inner_score_type == 'total':
                        components = data['components']
                        ddG = data['ddG']

                    elif score_type == 'noah_8,0A' and inner_score_type == 'positional':
                        ddG = data['ddG']
                    elif score_type == 'noah_8,0A' and inner_score_type == 'positional_twoscore':
                        ddG = data['ddG']
                    elif score_type == 'noah_8,0A' and inner_score_type == 'total':
                        ddG = data['ddG']
                    else:
                        continue
                        raise Exception('Unhandled score types: "%s", "%s".' % (score_type, inner_score_type))

                    ScoreMethodID = ScoreMethodMap[(score_type, inner_score_type)]
                    new_record = dict(
                        PredictionID = PredictionID,
                        ScoreMethodID = ScoreMethodID,
                        ScoreType = 'DDG',
                        StructureID = -1, # This score is for the Prediction rather than a structure
                        DDG = ddG,
                    )
                    assert(not(set(components.keys()).intersection(set(new_record.keys()))))
                    new_record.update(components)
                    ddGdb.insertDictIfNew('PredictionStructureScore', new_record, ['PredictionID', 'ScoreMethodID', 'ScoreType', 'StructureID'])

            if skip_if_we_have_pairs != None:
                # Skip this case if we have a certain number of existing records (much quicker since we do not have to extract the binary)
                num_wt = ddGdb.execute_select("SELECT COUNT(ID) AS NumRecords FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreType='WildType'", parameters=(PredictionID,))[0]['NumRecords']
                num_mut = ddGdb.execute_select("SELECT COUNT(ID) AS NumRecords FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreType='Mutant'", parameters=(PredictionID,))[0]['NumRecords']
                print(num_wt, num_mut)
                if num_wt == num_mut and num_mut == skip_if_we_have_pairs:
                    continue

            # Store the ddg_monomer scores for each structure
            grouped_scores = DDG_api.get_ddg_monomer_scores_per_structure(PredictionID)
            for structure_id, wt_scores in sorted(grouped_scores['WildType'].iteritems()):
                new_record = dict(
                    PredictionID = PredictionID,
                    ScoreMethodID = ScoreMethodMap[("kellogg", "total")],
                    ScoreType = 'WildType',
                    StructureID = structure_id,
                    DDG = None,
                )
                new_record.update(wt_scores)
                ddGdb.insertDictIfNew('PredictionStructureScore', new_record, ['PredictionID', 'ScoreMethodID', 'ScoreType', 'StructureID'])
            for structure_id, wt_scores in sorted(grouped_scores['Mutant'].iteritems()):
                new_record = dict(
                    PredictionID = PredictionID,
                    ScoreMethodID = ScoreMethodMap[("kellogg", "total")],
                    ScoreType = 'Mutant',
                    StructureID = structure_id,
                    DDG = None,
                )
                new_record.update(wt_scores)
                ddGdb.insertDictIfNew('PredictionStructureScore', new_record, ['PredictionID', 'ScoreMethodID', 'ScoreType', 'StructureID'])

            # Test to make sure that we can pick a best pair of structures (for generating a PyMOL session)
            assert(DDG_api.determine_best_pair(PredictionID) != None)



if __name__ == '__main__':
    DDG_api = ddGInterface()
    #ddGdb = ddgdbapi.ddGDatabase()
    determine_structure_scores(DDG_api)
