import sys
sys.path.insert(0, "..")
sys.path.insert(1, "../common")
import common.ddgproject
import common.colortext as colortext
from common.rosettahelper import ROSETTAWEB_SK_AAinv
from ddglib import help, dbapi
from ddglib.ddgfilters import *
from string import join
import pickle

ddGdb = common.ddgproject.ddGDatabase()

results = ddGdb.execute('''
SELECT ID, Chain, ResidueID, WildTypeAA, MutantAA, ddG, TIME_TO_SEC(TIMEDIFF(EndDate,StartDate))/60 as TimeTakenInMinutes FROM  `Prediction` 
INNER JOIN ExperimentMutation ON Prediction.ExperimentID = ExperimentMutation.ExperimentID
WHERE PredictionSet =  "lin-3K0NA"''')

colortext.message(len(results))

individual_results_by_position = {}
results_grouped_by_position = {}
wildtypes = {}
for r in results:
	assert(r["Chain"] == "A")
	assert(r["WildTypeAA"] != r["MutantAA"])
	resid = r["ResidueID"]
	ddG = pickle.loads(r["ddG"])["data"]["ddG"]
	
	individual_results_by_position[resid] = individual_results_by_position.get(resid) or {}
	individual_results_by_position[resid][r["MutantAA"]] = (r["WildTypeAA"], ddG, r["TimeTakenInMinutes"])
	
	wildtypes[resid] = r["WildTypeAA"]
	
	#results_grouped_by_position[resid] = results_grouped_by_position.get(resid) or []
	#results_grouped_by_position[resid].append((ddG, r["MutantAA"]))

F = open("PredictionsByMutation.csv", "w")
F.write("Chain\tResidueID\tWildType\tMutant\tddG\tWallTimeInMinutes\n")
# Warning: This sorting only works because there are no insertion codes (casting to int is okay)
for position in sorted(individual_results_by_position.keys(), key=lambda pos : int(pos)):
	for mutant, values in individual_results_by_position[position].iteritems():
		F.write("A\t%s\t%s\t%s\t%s\t%s\n" % (position, values[0], mutant, values[1], values[2]))
F.close()

aas = ROSETTAWEB_SK_AAinv.keys()
F = open("PredictionsByResidueID.csv", "w")
F.write("Chain\tResidueID\tWildType\tBestMutant\tBestMutant_ddG\t%s\n" % join(sorted(aas),"\t"))

# Warning: This sorting only works because there are no insertion codes (casting to int is okay)
minimum_best_mutant = 900
maximum_best_mutant = -900
best_mutants = {}
for position in sorted(individual_results_by_position.keys(), key=lambda pos : int(pos)):
	F.write("A\t%s\t%s\t" % (position, wildtypes[position]))
	
	results_grouped_by_position = [] 
	for mutant, values in individual_results_by_position[position].iteritems():
		results_grouped_by_position.append((values[1], mutant))
	
	best_mutant = sorted(results_grouped_by_position, key=lambda ddg_aa_pair : ddg_aa_pair[0])[0]
	best_mutants[position] = best_mutant[0] 
	minimum_best_mutant = min(minimum_best_mutant, best_mutant[0])
	maximum_best_mutant = max(maximum_best_mutant, best_mutant[0])
	F.write("%s\t%s\t" % (best_mutant[1], best_mutant[0]))
	
	results_grouped_by_position.append((0, wildtypes[position]))
	sorted_by_AA = sorted(results_grouped_by_position, key=lambda ddg_aa_pair : ddg_aa_pair[1])
	F.write(join([str(mtscore[0]) for mtscore in sorted_by_AA], "\t")) 
	F.write("\n")
F.close()

print('minimum_best_mutant', minimum_best_mutant)
print('maximum_best_mutant', maximum_best_mutant)

F = open("3K0NA_bfactors.pdb", "w")
pdbcontents = ddGdb.execute('''SELECT Content FROM Structure WHERE PDB_ID="3K0NA_lin"''')[0]["Content"].split("\n")
for line in pdbcontents:
	if line.startswith("ATOM  "):
		assert(line[21] == "A")
		assert(line[26] == " ")
		position = line[22:27].strip()
		newbfactor = "%.4f" % ((50.0 + (best_mutants[position] * 3.9))/100.0)
		assert(0 <= float(newbfactor) <= 1.0)
		assert(len(newbfactor) == 6)
		F.write("%s%s%s\n" % (line[0:60], newbfactor, line[66:]))
F.close()

	#ATOM    242  N   ALA A  33
