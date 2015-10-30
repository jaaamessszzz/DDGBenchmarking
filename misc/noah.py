# Script to get data for Noah - all double mutant data in the database and all single mutant data for those double mutants

import sys
import os
import json

sys.path.insert(0, "../..")
sys.path.insert(0, "..")
from ddglib import ddgdbapi
#from klab import colortext

def get_details(ExperimentID, expected_num_mutations = None):
    mutations = ddGdb.execute("SELECT * FROM ExperimentMutation WHERE ExperimentID=%s", parameters=(ExperimentID,))
    assert(expected_num_mutations == None or expected_num_mutations == len(mutations))
    mutation_list = []
    for m in mutations:
        d = {
            'Chain'         : m['Chain'],
            'WildTypeAA'    : m['WildTypeAA'],
            'ResidueID'     : m['ResidueID'],
            'MutantAA'      : m['MutantAA'],
        }
        mutation_list.append(d)

    DDG_data = list(ddGdb.execute("SELECT SecondaryID, Publication, Temperature, pH, Type, Value AS DDG FROM ExperimentAssayDDG INNER JOIN ExperimentAssay ON ExperimentAssay.ID=ExperimentAssayID WHERE ExperimentID=%s", parameters=(ExperimentID,)))
    return mutation_list, DDG_data

def get_double_mutants(ddGdb):
    double_mutant_ids = [(r['ID'], r['PDBFileID']) for r in ddGdb.execute("""
    SELECT ID, PDBFileID FROM (
        SELECT Experiment.ID, Experiment.PDBFileID, COUNT(ExperimentID) AS NumMutations
        FROM Experiment
        INNER JOIN(ExperimentMutation) ON ExperimentID=Experiment.ID
        GROUP BY ExperimentID
    ) AS t1 WHERE t1.NumMutations=2""")]

    json_object = {}
    for tpl in double_mutant_ids:
        ExperimentID = tpl[0]
        PDBFileID = tpl[1]
        mutation_list, DDG_data = get_details(ExperimentID, 2)
        json_object[ExperimentID] = {'PDB' : PDBFileID, 'Mutations' : mutation_list, 'DDGData' : DDG_data}
    return json_object

def get_list_of_single_mutants(ddGdb):
    single_mutant_ids = set([r['ID'] for r in ddGdb.execute("""
    SELECT ID, PDBFileID FROM (
        SELECT Experiment.ID, Experiment.PDBFileID, COUNT(ExperimentID) AS NumMutations
        FROM Experiment
        INNER JOIN(ExperimentMutation) ON ExperimentID=Experiment.ID
        GROUP BY ExperimentID
    ) AS t1 WHERE t1.NumMutations=1""")])
    return single_mutant_ids

def print_json(json_object):
    for ExperimentID, details in json_object.iteritems():
        print("> %d: %s" % (ExperimentID, details['PDB']))
        for m in details['Mutations']:
            print("$ %s %s %s %s" % (m['Chain'], m['WildTypeAA'], m['ResidueID'], m['MutantAA']))
        for ddg in details['DDGData']:
            pH = str(ddg['pH'] or '?')
            Temperature = str(ddg['Temperature'] or '?')
            print("@ %f, %s, %s, %s" % (ddg['DDG'], ddg['Type'], pH, Temperature))

def get_related_single_mutants(ExperimentID, details, single_mutant_ids, json_object):
    related_single_mutations = {}
    for m in details['Mutations']:
        related = ddGdb.execute("SELECT DISTINCT ExperimentID, PDBFileID FROM ExperimentMutation INNER JOIN Experiment ON Experiment.ID=ExperimentID WHERE Chain=%s AND ResidueID=%s AND WildTypeAA=%s AND MutantAA=%s", parameters=(m['Chain'], m['ResidueID'], m['WildTypeAA'], m['MutantAA']))
        for record in related:
            if record['ExperimentID'] in single_mutant_ids:
                RelatedExperimentID = record['ExperimentID']
                assert(RelatedExperimentID != ExperimentID)
                if record['PDBFileID'] == json_object[ExperimentID]['PDB']:
                    mutation_list, DDG_data = get_details(RelatedExperimentID, 1)
                    related_single_mutations[RelatedExperimentID] = {'PDB' : record['PDBFileID'], 'Mutations' : mutation_list, 'DDGData' : DDG_data}
    return related_single_mutations

if __name__ == '__main__':
    ddGdb = ddgdbapi.ddGDatabase()

json_object = get_double_mutants(ddGdb)
single_mutant_ids = get_list_of_single_mutants(ddGdb)

json_object_of_related_single_mutants = {}
for ExperimentID, details in json_object.iteritems():
    json_object_of_related_single_mutants[ExperimentID] = get_related_single_mutants(ExperimentID, details, single_mutant_ids, json_object)

count = {}
for ExperimentID, details in json_object_of_related_single_mutants.iteritems():
    num_related = len(details)
    count[num_related] = count.get(num_related, 0) + 1

count2 = {}
for ExperimentID, details in json_object.iteritems():
    count2[details['PDB']] = count2.get(details['PDB'], 0) + 1

for k, v in sorted(count2.iteritems(), key=lambda x:-x[1]):
    pass # print("%s: %d" % (k, v))



F = open('double_mutants.json', 'w')
F.write(json.dumps(json_object))
F.close()

F = open('related_single_mutants.json', 'w')
F.write(json.dumps(json_object_of_related_single_mutants))
F.close()
