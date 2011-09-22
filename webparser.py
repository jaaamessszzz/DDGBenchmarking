#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
from string import join
import math
import re
import db

#Rosetta3.1 release. mini revision 30964. minirosetta_database revision 30967.
#Rosetta3.2 release. mini revision 39284
#Rosetta3.3 release: mini revision 42942



fn = db.FieldNames()

def getPDB_IDs():
	#todo query db
	pass

def parseFloat(str):
	if str.strip() == "NaN":
		return None
	else:
		return float(str)

def parsePotapov():
	PredictionSet = "Potapov-2009"
	pdbIDs = getPDB_IDs()
	
	CCPBSA_ID = None 
	EGAD_ID =  None
	FoldX_ID =  None
	Hunter_ID = None
	IMutant2_ID = None
	Rosetta_ID =  None 
	
	#todo: Get indices for tools above
	#todo: Only add data if none is in db for Potapov-2009

	Prediction = {fn.PredictionSet	: PredictionSet}
	
	mutantsfile = os.path.join("rawdata", "mutants.txt")
	F = open(mutantsfile)
	for line in F.read().split("\n"):
		if line.strip() and line[0] != '#':
			
			Mutations = []
			ExperimentChains = []
			ExperimentalScores = []
			chainID = line[9]
			pdbid = line[5:9]
			chain = line[25] # IgnoredField 
			
			Structure = {fn.PDB_ID : pdbid} 
			#insertStructure(Structure, pdbIDs) # todo
			
			Experiment = {
				fn.Structure	: pdbid,
				fn.Mutant		: None,
				fn.Source		: "Potapov",
				fn.SourceID		: int(line[0:4])
			}	
			
			Experiment["Mutations"] = [{
				fn.Chain 		: chainID,
				fn.WildTypeAA	: line[14],
				fn.MutantAA		: line[19],
				fn.ResidueID	: int(line[20:25])
				}]
			
			Experiment["ExperimentChains"] = [{
				fn.Chain		: chainID
				}]
			
			Experiment["ExperimentalScores"] = [{
				fn.ddG				: parseFloat(line[26:35]),
				fn.NumberOfMeasurements : int(line[36:41])
				}]
			
			ExperimentID = None
			#ExperimentID = insertExperiment(Experiment) #todo
			Prediction[fn.ExperimentID] = ExperimentID
			
			Prediction[fn.ToolID] = CCPBSA_ID
			Prediction[fn.ddG] = parseFloat(line[41:51]) 
			Prediction[fn.Description] = {"MutationUsedInEvaluatingTheMethod" : int(line[51:53]) == 1}
			#insertPrediction(Prediction) #todo
			
			Prediction[fn.ToolID] = EGAD_ID
			Prediction[fn.ddG] = parseFloat(line[53:63]) 
			Prediction[fn.Description] = {"MutationUsedInEvaluatingTheMethod" : int(line[63:65]) == 1}
			#insertPrediction(Prediction) #todo 
			
			Prediction[fn.ToolID] = FoldX_ID
			Prediction[fn.ddG] = parseFloat(line[65:75]) 
			Prediction[fn.Description] = {"MutationUsedInEvaluatingTheMethod" : int(line[75:77]) == 1}
			#insertPrediction(Prediction) #todo

			Prediction[fn.ToolID] = Hunter_ID
			Prediction[fn.ddG] = parseFloat(line[77:87]) 
			Prediction[fn.Description] = {"MutationUsedInEvaluatingTheMethod" : int(line[87:89]) == 1}
			#insertPrediction(Prediction) #todo

			Prediction[fn.ToolID] = IMutant2_ID
			Prediction[fn.ddG] = parseFloat(line[89:99]) 
			Prediction[fn.Description] = {"MutationUsedInEvaluatingTheMethod" : int(line[99:101]) == 1}
			#insertPrediction(Prediction) #todo

			Prediction[fn.ToolID] = Rosetta_ID
			Prediction[fn.ddG] = parseFloat(line[101:111]) 
			Prediction[fn.Description] = {"MutationUsedInEvaluatingTheMethod" : int(line[111:113]) == 1}
			#insertPrediction(Prediction) #todo
					
	F.close()

def parseProTherm():
	pdbIDs = getPDB_IDs()
	
	ID = None
	mutation = {}
	
	listOfNoIDs = []
	listOfNoLengths = []
	mutantProcessingErrors = []
	ddGProcessingErrors = []
	chainProcessingErrors = []
	
	experiments = {}
	totalcount = 0
	count = 0
	MAX_NUMRES = 350
	MAX_STANDARD_DEVIATION = 0.5
	chains = {}
	singleErrors = 0
	
	# These fields of ProTherm records cannot be empty for our purposes
	requiredFields = ["NO.", "PDB_wild", "LENGTH", "ddG", "MUTATION", "MUTATED_CHAIN"]
	
	newlist = []
	protherm = os.path.join("rawdata", "ProTherm.dat")
	F = open(protherm)
	while(True):
		# Read a record
		record = {}
		line = F.readline()
		while line and not(line.startswith("//")):
			if line[0] != "*" and line[0] != " ":
				#record.append(line)
				line = line.split()
				if len(line) > 1:
					record[line[0]] = line[1:]
				else:
				 	record[line[0]] = None
			line = F.readline()
		
		# Parse the results
		if record:
			totalcount += 1	
			
			# Find out whether we have enough information
			store = True
			ID = int(record["NO."][0])
			missingFields = []
			for field in requiredFields:
				if not record[field]:
					store = False
					missingFields.append(field)
			if not store:
				if len(missingFields) == 1:
					if record["MUTATION"][0] and record["MUTATION"][0] != "wild" and missingFields[0] != "ddG":
						if not record["MUTATED_CHAIN"]:
							print("Error processing chain: ID %d, no chain" %  (ID))
						elif not record["LENGTH"]:
							print("Error processing length: ID %d, no length" %  (ID))
						elif not record["PDB_wild"]:
							print("Error processing PDB ID: ID %d, no PDB ID" %  (ID))
						else:
							print("Error processing structure: ID %d, no %s " % (ID, missingFields[0]))
						singleErrors += 1
				continue
			else:
				# Only allow proteins with <= MAX_NUMRES residues
				if int(record["LENGTH"][0]) > MAX_NUMRES:
					continue
			
			pdbid = record["PDB_wild"][0]
			
			# Parse chain
			if len(record["MUTATED_CHAIN"]) == 1:
				chainID = record["MUTATED_CHAIN"][0]
				chains[chainID] = True
			else:
				print("Error processing chain: ID %d, %s" %  (ID, record["MUTATED_CHAIN"]))
				store = False
			
			# Parse mutation
			mutationline = record["MUTATION"]
			if len(mutationline) == 3:
				try:
					ResidueID = int(mutationline[1])
					WildTypeAA = mutationline[0]
					MutantAA = mutationline[2]
				except:
					print("Error processing mutation: ID %d, %s" % (ID, record["MUTATION"])) 
					store = False
			else:
				store = False
				
			# Parse ddG
			try:
				ddGline = join(record["ddG"], " ")
				ddG = float(ddGline)				
			except:
				idx = ddGline.find("kJ/mol")
				if idx != -1:
					try:
						ddG = float(ddGline[0:idx])
					except:
						print("Error processing ddG: ID %d, %s" % (ID, record["ddG"]))
						store = False
				else:
					idx = ddGline.find("kcal/mol")
					if idx != -1:
						try:
							ddG = 4.184 * float(ddGline[0:idx])
						except:
							print("Error processing ddG: ID %d, %s" % (ID, record["ddG"]))
							store = False
					else:
						print("Error processing ddG: ID %d, %s" % (ID, record["ddG"]))
						store = False 
				
				
			if store:
				# Parse mutant
				mutant = None
				if record["PDB_mutant"]:
					mutant = record["PDB_mutant"][0]
				
				Structure = {fn.PDB_ID : pdbid} 
				#insertStructure(Structure, pdbIDs) # todo
				
				# NOTE: THIS IS A PATCH FOR BOTH BAD AND MISMATCHED DATA IN ProTherm 
				if pdbid in ["1CSP", "2LZM"]:
					chainID = "A"
					
				Experiment = {
					fn.Structure	: pdbid,
					fn.Mutant		: mutant,
					fn.Source		: "ProTherm-2008-09-08-23581",
					fn.SourceID		: ID
				}	
				
				Experiment["Mutations"] = [{
					fn.Chain 		: chainID,
					fn.WildTypeAA	: WildTypeAA,
					fn.MutantAA		: MutantAA,
					fn.ResidueID	: ResidueID
					}]
				
				Experiment["ExperimentChains"] = [{
					fn.Chain		: chainID
					}]
				
				Experiment["ExperimentalScores"] = [{
					fn.ddG					: ddG,
					fn.NumberOfMeasurements : 1
					}]
				
									
				# All necessary information is present. Store the record.
				mutationID = "%s-%s-%s%d%s" % (pdbid, chainID, WildTypeAA, ResidueID, MutantAA)
				if experiments.get(mutationID):
					for e in experiments[mutationID]:
						set1 = set([k[fn.Chain] for k in e["ExperimentChains"]])
						set2 = set([k[fn.Chain] for k in Experiment["ExperimentChains"]])
						if set1.difference(set2):
							
							print(record["LENGTH"][0])
							print(Experiment)
							print(e)
							print(pdbid)
				experiments[mutationID] = experiments.get(mutationID) or [] 
				experiments[mutationID].append(Experiment)
				count += 1
				newlist.append(ID)
		else:
			break	

	# Exclude experiments with high standard deviation #todo	
	highvariancecount = 0
	for mutation, experimentalResults in experiments.iteritems():
		n = len(experimentalResults)
		if n > 1:
			n = float(n)
			sum = 0
			for experiment in experimentalResults:
				sum += experiment["ExperimentalScores"][0]["ddG"]
			mean = sum / n
			squaredsum = 0
			for experiment in experimentalResults:
				diff = (experiment["ExperimentalScores"][0]["ddG"] - mean)
				squaredsum += diff * diff
			stddev = math.sqrt(squaredsum / n)
			if stddev > MAX_STANDARD_DEVIATION:
				highvariancecount += 1
				#print(mutation)
				#print("\tNumber of results: %d" % n) #result["ddG"])
				#print("\tStandard deviation: %0.2f" % stddev) #result["ddG"])
			else:
				Experiment = experimentalResults[0]
				# todo: We should sanity-check here and make sure all common fields are equal 
				if False:
					Experiment = {
						fn.Structure	: pdbid,
						fn.Mutant		: None,
						fn.Source		: "Potapov",
						fn.SourceID		: int(line[0:4])
					}	
					
					Experiment["Mutations"] = [{
						fn.Chain 		: chainID,
						fn.WildTypeAA	: line[14],
						fn.MutantAA		: line[19],
						fn.ResidueID	: int(line[20:25])
						}]
					
					Experiment["ExperimentChains"] = [{
						fn.Chain		: chainID
						}]
				
				Experiment["ExperimentalScores"] = [{
					fn.ddG					: mean, # Store the mean of all measurements
					fn.NumberOfMeasurements : int(n)
					}]
				#ExperimentID = insertExperiment(Experiment) #todo
			
		elif n == 1:
			Experiment = experimentalResults[0]
			#ExperimentID = insertExperiment(Experiment) #todo
			
	
	print("")
	print("Number of mutations: %d" % totalcount)
	print("Number of acceptable mutations: %d" % count)
	print("Number of unique mutations: %d" % len(experiments))
	print("Number of potentially correctable single errors: %d" % singleErrors)
	
	PDBIDs = [id[0:4] for id in experiments.keys()]
	pd = {}
	for x in PDBIDs:
		pd[x] = True
	print("Number of unique PDB IDs: %d" % len(pd))
	print("Chain IDs: %s" % join(sorted(chains.keys()), ", "))
	
	#print("listOfNoIDs",listOfNoIDs)
	errors = mutantProcessingErrors + ddGProcessingErrors + chainProcessingErrors 
	for e in errors:
		print(e)
	
	for k,v in experiments.iteritems():
		#print(k)
		for i in range(len(v)):
			pass
			#print(i, v[i])
	
	print("Number of mutations with standard deviation > %0.2f: %d" % (MAX_STANDARD_DEVIATION, highvariancecount))

	
	
	print("***\n\n")
	#print("listOfNoLengths",listOfNoLengths)
	F = open("PDBIDs.txt", "w")
	F.write(join(PDBIDs, "\n"))
	F.close()
	

def parseUniProt():
	uniprot = os.path.join("rawdata", "uniprotmapping.txt")
	F = open(uniprot)
	lines = F.read().split("\n")[1:]
	F.close()
	
	mapping = {}
	for line in lines:
		PDBID, ID = line.split("\t")
		#print(PDBID, ID)
		if mapping.get(ID):
			print("%s -> %s, %s -> %s" % ( mapping[ID], ID, PDBID, ID))
		mapping[ID] = PDBID
	print(len(mapping.keys()))

def parseRawData():
	parsePotapov()
	parseProTherm()
	parseUniProt()
	
def main():
	parseRawData()
	sys.exit(0)
	#All results were produced with revision 32231 of rosetta, and revision 32257 of the rosetta database. 
	commonstr = "in:file:s <INPUT_PDB> -resfile <RESFILE> -database <DATABASE> -ignore_unrecognized_res –in:file:fullatom –constraints::cst_file <CONSTRAINTS_FILE>"
	softrep = " -score:weights soft_rep_design"
	hardrep = "-score:weights standard –score:patch score12"
	minnohardrep = " -ddg::minimization_scorefunction standard –ddg::minimization_patch score12"
	
	executable = "fix_bb_monomer_ddg.linuxgccrelease"
	
	protocols = range(0,21)
	Description = "Fixed backbone, sidechain repacking extent: 1 residue"
	protocols[1] = [#{
		#"Description" : "Fixed backbone",
		#"CommandLine" : [
			executable, 
			'-ddg::weight_file', 'soft_rep_design',
			'-ddg::iterations', '1', 
			'-ddg::local_opt_only', 'true',
			'-ddg::min_cst', 'false', 
			'-ddg::mean', 'false',
			'-ddg::min', 'true', 
			'-ddg::sc_min_only', 'false',
			'-ddg::opt_radius', '0.1']#}
	
	#protocols[2] = [
	
	
	
	Description = "Fixed backbone"
	protocols[3] = [executable, 
					'-ddg::weight_file', 'soft_rep_design',
					'-ddg::iterations', '50',
					'-ddg::local_opt_only', 'true',
					'-ddg::min_cst', 'false',
					'-ddg::mean', 'true',
					'-ddg::min', 'false',
					'-ddg::sc_min_only', 'false',
					'-ddg::opt_radius', '8.0',
					'-ddg::ramp_repulsive', 'false'] 
	protocols[5] = [executable, 
					'-ddg::weight_file', 'soft_rep_design',
					'-ddg::iterations', '50',
					'-ddg::local_opt_only', 'true',
					'-ddg::min_cst', 'true',
					'-ddg::mean', 'false',
					'-ddg::min', 'true',
					'-ddg::sc_min_only', 'true',
					'-ddg::opt_radius', '8.0'
					'-ddg::ramp_repulsive', 'false', 
					'-ddg::minimization_scorefunction', 'standard',
					'-ddg::minimization_patch', 'score12']
					
	protocols[6] = [executable,
					'-ddg::weight_file', 'soft_rep_design',
					'-ddg::iterations', '50', 
					'-ddg::local_opt_only', 'false',
					'-ddg::min_cst', 'false',
					'-ddg::mean', 'true',
					'-ddg::min', 'false',
					'-ddg::sc_min_only', 'false']
	
	protocols[8] = [executable, 
					'-ddg::weight_file', 'soft_rep_design',
					'-ddg::iterations', '50',
					'-ddg::local_opt_only', 'false',
					'-ddg::min_cst', 'true',
					'-ddg::mean', 'false', 
					'-ddg::min', 'true', 
					'-ddg::sc_min_only', 'true', 
					'-ddg::ramp_repulsive', 'false', 
					'-ddg::minimization_scorefunction', 'standard',
					'-ddg::minimization_patch', 'score12']
	###
	protocols[10] = [executable, 
					'-ddg::weight_file', 'standard_plus_score12.wts',
					'-ddg::iterations', '1',
					'-ddg::local_opt_only', 'true',
					'-ddg::min_cst', 'true',
					'-ddg::mean', 'false',
					'-ddg::min', 'true',
					'-ddg::sc_min_only', 'false',
					'-ddg::opt_radius', '0.1',
					'-ddg::ramp_repulsive', 'true',
					'-ddg::minimization_scorefunction', 'standard',
					'-ddg::minimization_patch', 'score12']
	protocols[13] = [executable, 
					'-ddg::weight_file', 'soft_rep_design',
					'-ddg::iterations', '50',
					'-ddg::local_opt_only', 'true',
					'-ddg::min_cst', 'true',
					'-ddg::mean', 'false',
					'-ddg::min', 'true',
					'-ddg::sc_min_only', 'false',
					'-ddg::opt_radius', '8.0',
					'-ddg::ramp_repulsive', 'true', 
					'-ddg::minimization_scorefunction', 'standard',
					'-ddg::minimization_patch', 'score12']
	protocols[16] = [executable, 
					'-ddg::weight_file', 'soft_rep_design',
					'-ddg::iterations', '50',
					'-ddg::local_opt_only', 'false',
					'-ddg::min_cst', 'true',
					'-ddg::mean', 'false',
					'-ddg::min', 'true',
					'-ddg::sc_min_only', 'false',
					'-ddg::ramp_repulsive', 'true', 
					'-ddg::minimization_scorefunction', 'standard',
					'-ddg::minimization_patch', 'score12']
	protocols[18] = protocols[16]
	protocolsconstraints18 = ['make_cst_file.linuxiccrelease',
							'\xe2\x80\x93ddg::distance_from_mutsite', '10.0',
							'\xe2\x80\x93ddg::strict_cst', '0.5',
							'-ddg::loose_cst', '2']
	protocols[19] = protocols[16]
	protocolsconstraints19 = None
	###
	protocols[20] = ['ensemble_generator_score12_sidechain_ver2.linuxgccrelease',
					'-sc_min_only', 'false', 
					'-ddg::ramp_repulsive', 'true',
					'-ddg::constraint_weight', '1.0', 
					'-nstruct', '200', 
					'-ddg::min_with_cst', 'true',
					'-ddg::temperature', '10.0',
					'-ddg::use_bound_cst', 'true']
	
	i = 0
	for p in protocols:
		if type(p) == type(""):
			print(i, p.split(" "))
		i += 1	
	#print(protocolsconstraints18.split(" "))
	#'protocols[1] = "
	
	
main() 
