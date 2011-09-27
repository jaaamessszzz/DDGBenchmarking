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

ExperimentSet = db.ExperimentSet
Prediction = db.Prediction
PDBStructure = db.PDBStructure
fn = db.FieldNames()

AllPDBIDs = {}

def downloadPDBs(pdbIDlist):
	for pdbID in pdbIDlist:
		db.downloadPDB(pdbID)
	pass

def parseFloat(str):
	if str.strip() == "NaN":
		return None
	else:
		return float(str)

def pdb2toRes(pdbID):
	''' Used for debugging.'''
	ROSETTAWEB_SK_AA = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
					"HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
					"PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
					"TRP": "W", "TYR": "Y"}
	mutantsfile = os.path.join("pdbs", "%s.pdb" % pdbID)
	F = open(mutantsfile)
	
	lastres = None
	resnumber = 1
	
	for line in F.read().split("\n"):
		if line[0:6] == "ATOM  ":
			if line[22:27] != lastres: # include insertion code
				print("%s\t%s\t%s" % (line[22:27].strip(), line[17:20], line[21]))
				resnumber += 1
				lastres = line[22:27] 
	F.close()
	
def pdb2sequentialAtoms(pdbID, chainsToParse = []):
	ROSETTAWEB_SK_AA = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
					"HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
					"PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
					"TRP": "W", "TYR": "Y"}
	mutantsfile = os.path.join("pdbs", "%s.pdb" % pdbID)
	F = open(mutantsfile)
	
	lastres = None
	resnumber = 1
	
	resmap = {}
	for line in F.read().split("\n"):
		if line[0:6] == "ATOM  ":
			if not chainsToParse or line[21] in chainsToParse:
				if line[22:27] != lastres: # include insertion code
					resmap[resnumber] = (line[22:27].strip(), ROSETTAWEB_SK_AA[line[17:20]], line[21])
					resnumber += 1
					lastres = line[22:27] 
	F.close()
	return resmap
	
def convertSensDatabase(ddGdb):
	'''This function is used to convert Sen's spreadsheet to a uniform spreadsheet.'''
	mutantsfile = os.path.join("rawdata", "sens-complex-dataset-exp-ddg.csv")
	outputfile = os.path.join("rawdata",  "sens-complex-dataset-exp-ddg-reformatted.csv")
	headers = ['Pair', 'Wtpdb', 'MUTpdb', 'wt', 'Sites', 'mut', 'Tmp', 'ddGexp']
	recordlength = len(headers)
	newheaders = ['Interface', 'Structure', 'Mutant' ,'Chains', 'Mutations', 'ddGexp']
	#print("Citations")
	#print("BPTI-BCHYM\thttp://dx.doi.org/10.1002/pro.5560060902")

	missingddGs = []
	
	#0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#QVQLQQSGTELVKSGASVKLSCTASGFNIKDTHMNWVKQRPEQGLEWIGRIDPANGNIQYDPKFRGKATITADTSSNTAYLQLSSLTSEDTAVYYCATKVIYYQGRGAMDYWGQGTTLTVS
	mapping = {
		"BPTI-BCHYM" : ("1CBW", ['D']),
		"CHTP-OMTKY3" : ("1CHO", ['G']),
		"D1.3-E5.2" : ("1DVF", ['A','B','C','D']),
		"E9-IM9" : ("1EMV", ['A']), # B pruned unnecessarily
		"FV-HEL" : ("2DQJ", ['L', 'H', 'Y']),
		"IGG1-HEWL" : ("1VFB", ['A', 'B', 'C']),
		"IL4-IL4Ra" : ("1IAR", ['A','B']),
		"Rnase-Ang" : ("1A4Y", ['B','D']), # A pruned necessarily, E pruned unnecessarily
		"TEM1-BLIP" : ("1JTG", ['A','B']) # C and D pruned unnecessarily
		}
	pdbInfo = {}
	
	for interface, data in mapping.iteritems():
		pdbID = data[0]
		chainsToParse = data[1]
		pdbInfo[interface] = pdb2sequentialAtoms(pdbID, chainsToParse)
		
	D13E52regex = re.compile("(\d+)(\w)(\d+)(\w)")
	FVHELregex = re.compile("(\d+\w)")
	TEM1BLIPregex1 = re.compile("^\d+_\w$")
	TEM1BLIPregexn = FVHELregex
	
	F = open(mutantsfile)
	contents = F.read().split("\n")
	F.close()
	
	outfile = open(outputfile, "w")
	#print(newheaders)
	outfile.write(join(newheaders, "\t") + "\n")
	for i in range(1, len(contents)):
		line = contents[i].split("\t")
		lineno = i+1
		length = len(line)
		assert((length == 1 and not(line[0])) or (length == recordlength))
		if length == recordlength:
			try:
				ddG = float(line[7])
			except:
				missingddGs.append("Missing ddG on line %d.\n\t%s" % (lineno, line))
				ddG = None
				continue
			if line[2].strip() == "--":
				line[2] = None
			mutant = line[2] or ""
			interface = line[0].strip()
			pdbID = mapping.get(interface,["",""])[0]
			resmap = pdbInfo.get(interface)
			
			if interface == "BPTI-BCHYM":
				chains = "D"
				chain = 'D'
				wt = line[3]
				resid = int(line[4])
				mut = line[6].split("_")[1]
				
				#Sanity check
				assert(resmap[resid][1] == wt)
				resid = resmap[resid][0]
				
				mutations = '%s-%s%s%s' % (chain, wt, resid, mut)
				newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
				#print(newline) 
				outfile.write(newline + "\n")
			elif interface == "CHTP-OMTKY3":
				chains = "G"
				chain = 'G'
				wt = line[3].strip()
				resid = int(line[4])
				mut = line[5].strip()
				
				#Sanity check
				assert(resmap[resid][1] == wt)
				resid = resmap[resid][0]
				
				mutations = '%s-%s%s%s' % (chain, wt, resid, mut)
				newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
				#print(newline) 
				outfile.write(newline + "\n")
			elif interface == "D1.3-E5.2":
				chains = "A,B,C,D"
				wts = line[3]
				if len(wts) == 3:
					result = D13E52regex.match(line[6])
					if not result:
						raise Exception("Regex for D1.3-E5.2 failed.\n%s" % line)
					muts = [result.group(2), result.group(4)]
					resids = [int(result.group(1)), int(result.group(3))]
					
					wts = wts.strip().split(",")
					if lineno == 72:
						wts[1] = "H"
					elif lineno == 73:
						wts[1] = "D"
					elif lineno == 63 or lineno == 68 or lineno == 70:
						resids[1] += 1 
					assert(len(wts) == 2)
					
					#Sanity check
					try:
						assert(resmap[resids[0]][1] == wts[0])
					except:
						print("\nMismatch on line %d" % lineno)
						print(line)
						print(resids[0])
						print("Read wildtype: %s" % resmap[resids[0]][1])
						print("Sen's wildtype: %s\n" % wts[0])
						#raise Exception()
					try:
						assert(resmap[resids[1]][1] == wts[1])
					except:
						print("\nMismatch on line %d" % lineno)
						print(line)
						print(resids[1])
						print("Read wildtype: %s" % resmap[resids[1]][1])
						print("Sen's wildtype: %s\n" % wts[1])
						#raise Exception()
					chain  = [resmap[resids[0]][2], resmap[resids[1]][2]]
					resids = [resmap[resids[0]][0], resmap[resids[1]][0]]
					
					mutations = ['%s-%s%s%s' % (chain[k], wts[k], resids[k], muts[k]) for k in [0,1]]
					newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, join(mutations,","), ddG)
					#print(newline) 
					outfile.write(newline + "\n")
				elif len(wts) == 1:
					mutchain = line[6].split("_")
					resid = int(mutchain[0])
					mut = mutchain[1]
					
					if lineno == 81 or lineno == 82:
						resid += 1 
					
					#Sanity check
					try:
						assert(resmap[resid][1] == wts)
					except:
						print("\nMismatch on line %d" % lineno)
						print(line)
						print(resid)
						print("Read wildtype: %s" % resmap[resid][1])
						print("Sen's wildtype: %s\n" % wts)
						
					chain  = resmap[resid][2]
					resid = resmap[resid][0]

					mutations = '%s-%s%s%s' % (chain, wts, resid, mut)
					newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
					#print(newline) 
					outfile.write(newline + "\n")				
				else:
					raise Exception("Unhandled case in D1.3-E5.2.\n%s" % line)
			elif interface == "E9-IM9":
				chains = "A"
				chain = 'A'
				wt = line[3].strip()
				resid = int(line[4])
				mut = line[5].strip()
				
				#Sanity check
				assert(resmap[resid][1] == wt)
				resid = resmap[resid][0]
				
				mutations = '%s-%s%s%s' % (chain, wt, resid, mut)
				newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
				#print(newline) 
				outfile.write(newline + "\n")
			elif interface == "FV-HEL" or interface == "IGG1-HEWL":
				if interface == "FV-HEL":
					chains = "L,H,Y"
				elif interface == "IGG1-HEWL":
					chains = "A,B,C"
				wt = line[3].strip()
				if len(wt) == 1:
					resid = int(line[4])
					mut = line[5].strip()
					
					#Sanity check
					assert(resmap[resid][1] == wt)
					chain = resmap[resid][2]
					resid = resmap[resid][0]
					
					mutations = '%s-%s%s%s' % (chain, wt, resid, mut)
					newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
					#print(newline) 
					outfile.write(newline + "\n")
				else:
					wts = wt.split(",")
					
					result = FVHELregex.findall(line[6])
					if not result:
						raise Exception("Regex for D1.3-E5.2 failed.\n%s" % line)
					assert(len(wts) == len(result))
					
					# Dirty hack to get around typos
					if interface == "IGG1-HEWL":
						if lineno == 142:
							assert(line[4] == "32347")
							result =["32A", "347A"]
						elif lineno == 143:
							assert(line[4] == "32344")
							result =["32A", "344A"]
						elif lineno == 151:
							assert(line[4] == "92347")
							result =["92A", "347A"]
						elif lineno == 152:
							assert(line[4] == "92352")
							result =["92A", "352A"]
						elif lineno == 153:
							assert(line[4] == "92344")
							result =["92A", "344A"]
						elif lineno == 154:
							assert(line[4] == "92348")
							result =["92A", "348A"]
					
					#Sanity check
					mutations = []
					for k in range(len(wts)):
						wt = wts[k]
						resid = int(result[k][0:-1])
						mut = result[k][-1]
						
						try:
							assert(resmap[resid][1] == wt)
						except:
							print("\nMismatch on line %d" % lineno)
							print(line)
							print(resid)
							print("Read wildtype: %s" % resmap[resid][1])
							print("Sen's wildtype: %s\n" % wt)
						
						
						chain = resmap[resid][2]
						resid = resmap[resid][0]
						mutations.append('%s-%s%s%s' % (chain, wt, resid, mut))
					mutations = join(mutations, ",")
					newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
					#print(newline) 
					outfile.write(newline + "\n")
			elif interface == "IL4-IL4Ra":
				chains = "A,B"
				
				wt = line[3]
				resid = int(line[4])
				mut = line[6].split("_")[1]
				
				#Sanity check
				assert(resmap[resid][1] == wt)
				chain = resmap[resid][2]
				resid = resmap[resid][0]
				
				mutations = '%s-%s%s%s' % (chain, wt, resid, mut)
				newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
				#print(newline) 
				outfile.write(newline + "\n")
			elif interface =="Rnase-Ang":
				chains = "B,D"
				wt = line[3].strip()
				if len(wt) == 1:
					resid = int(line[4][1:-1])
					mut = line[5].strip()
					
					#Sanity check
					assert(resmap[resid][1] == wt)
					chain = resmap[resid][2]
					resid = resmap[resid][0]
					
					mutations = '%s-%s%s%s' % (chain, wt, resid, mut)
					newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
					#print(newline) 
					outfile.write(newline + "\n")
				else:
					wts = wt.split(",")
					
					result = FVHELregex.findall(line[6])
					if not result:
						raise Exception("Regex for D1.3-E5.2 failed.\n%s" % line)
					assert(len(wts) == len(result))
					
					# Dirty hack to get around typos
					if lineno == 239:
						assert(line[4] == "R5AD558A")
						result =["5A", "558A"]
					elif lineno == 240:
						assert(line[4] == "R5AY557A")
						result =["5A", "557A"]
						
					#Sanity check
					mutations = []
					for k in range(len(wts)):
						wt = wts[k]
						resid = int(result[k][0:-1])
						mut = result[k][-1]
						
						try:
							assert(resmap[resid][1] == wt)
						except:
							print("\nMismatch on line %d" % lineno)
							print(line)
							print(resid)
							print("Read wildtype: %s" % resmap[resid][1])
							print("Sen's wildtype: %s\n" % wt)
						
						
						chain = resmap[resid][2]
						resid = resmap[resid][0]
						mutations.append('%s-%s%s%s' % (chain, wt, resid, mut))
					mutations = join(mutations, ",")
					newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
					#print(newline) 
					outfile.write(newline + "\n")
			elif interface =="TEM1-BLIP":
				chains = "A,B"
				wt = line[3].strip()
				if not wt:
					wt = None
				onemut = TEM1BLIPregex1.match(line[6])
				if onemut:
					resmut = line[6].split("_")
					resid = int(resmut[0])
					mut = resmut[1]
					
					# Awful hack to get around what seems to be a counting error in Sen's spreadsheet
					if resid >= 189:
						resid -= 1
						
					#Sanity check
					if wt:
						try:
							assert(resmap[resid][1] == wt)
						except:
							print("\nMismatch on line %d" % lineno)
							print(line)
							print(resid)
							print("Read wildtype: %s" % resmap[resid][1])
							print("Sen's wildtype: %s\n" % wt)
							
					wt = resmap[resid][1]
					chain = resmap[resid][2]
					resid = resmap[resid][0]
					
					mutations = '%s-%s%s%s' % (chain, wt, resid, mut)
					newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
					#print(newline) 
					outfile.write(newline + "\n")
				else:
					result = TEM1BLIPregexn.findall(line[6])
					if not result:
						raise Exception("Regex for D1.3-E5.2 failed.\n%s" % line)
					
					#Sanity check
					mutations = []
					assert(not(wt))
					for k in range(len(result)):
						resid = int(result[k][0:-1])
						
						# Awful hack to get around what seems to be a counting error in Sen's spreadsheet
						if resid >= 189:
							resid -= 1
						
						wt = resmap[resid][1]
						mut = result[k][-1]
						
						chain = resmap[resid][2]
						resid = resmap[resid][0]
						mutations.append('%s-%s%s%s' % (chain, wt, resid, mut))
					mutations = join(mutations, ",")
					newline = "%s\t%s\t%s\t%s\t%s\t%f" % (interface, pdbID, mutant, chains, mutations, ddG)
					#print(newline) 
					outfile.write(newline + "\n")				
	
	outfile.close()
			
	if missingddGs:
		print("\n**Errors parsing Sen's dataset")
		for e in missingddGs:
			print(e)
	
	# Sanity check
	infile = open(mutantsfile, "r")
	outfile = open(outputfile, "r")
	incontents = infile.read().split("\n")
	outcontents = outfile.read().split("\n")
	infile.close
	outfile.close
	assert(len(incontents) == len(outcontents) + len(missingddGs))
	incontents.pop(204) # Hack for missing ddG
	for i in range(1, len(incontents)):
		if incontents[i]:
			assert(outcontents[i])
			inline = incontents[i].split("\t")
			outline = outcontents[i].split("\t")
			if inline:
				assert(float(inline[7]) == float(outline[5]))
		else:
			assert(not(outcontents[i]))
			
				
	

def parseSensDatabase(ddGdb):
	PredictionSet = "SenLiuComplexDatasetExp"
	mutantsfile = os.path.join("rawdata", "sens-complex-dataset-exp-ddg-reformatted.csv")
	headers = ['Interface', 'Structure', 'Mutant' ,'Chains', 'Mutations', 'ddGexp']
	recordlength = len(headers)
	F = open(mutantsfile)
	linecount = 1
	for line in F.read().split("\n")[1:]:
		line = line.split("\t")
		length = len(line)
		assert((length == 1 and not(line[0])) or (length == recordlength))
		if length == recordlength:
			pdbID = line[1].upper()
			mutant = line[2]
			AllPDBIDs[pdbID] = True
			if mutant:
				AllPDBIDs[mutant] = True
				
			interface = line[0]
			chainIDs = line[3].split(",")
			mutations = line[4].split(",")
			ddG = float(line[5])
			Structure = PDBStructure(pdbID)
			Structure.commit(ddGdb)
			
			Experiment = ExperimentSet(pdbID, "SenLiu", interface = interface)
			if mutant:
				Experiment.addMutant(mutant)
			for mutation in mutations:
				chainID = mutation[0]
				wildtypeAA = mutation[2]
				residueID = mutation[3:-1]
				mutantAA = mutation[-1]
				Experiment.addMutation(chainID, residueID, wildtypeAA, mutantAA)
			for chainID in chainIDs:
				 Experiment.addChain(chainID)
			Experiment.addExperimentalScore(linecount, ddG)
			Experiment.commit(ddGdb)
			
			ExperimentID = 0 #None#Experiment.getDBID() todo
			linecount += 1
	F.close()
	
def parsePotapov(ddGdb):
	PredictionSet = "Potapov-2009"
	
	CCPBSA_ID = None or 1 # todo
	EGAD_ID =  None or 1 # todo
	FoldX_ID =  None or 1 # todo
	Hunter_ID = None or 1 # todo
	IMutant2_ID = None or 1 # todo
	Rosetta_ID =  None  or 1 # todo
	
	#todo: Get indices for tools above
	#todo: Only add data if none is in db for Potapov-2009

	mutantsfile = os.path.join("rawdata", "mutants.txt")
	F = open(mutantsfile)
	for line in F.read().split("\n"):
		if line.strip() and line[0] != '#':
			
			chainID = line[9]
			pdbID = line[5:9].upper()
			AllPDBIDs[pdbID] = True
			
			Structure = PDBStructure(pdbID)
			Structure.commit(ddGdb)
			
			Experiment = ExperimentSet(pdbID, "Potapov")
			Experiment.addMutation(chainID, int(line[20:25]), line[14], line[19])
			Experiment.addChain(chainID)
			Experiment.addExperimentalScore(int(line[0:4]), parseFloat(line[26:35]), int(line[36:41]))
			Experiment.commit(ddGdb)
			
			ExperimentID = 0 #None#Experiment.getDBID() todo
			
			p = Prediction(ExperimentID, PredictionSet, CCPBSA_ID, parseFloat(line[41:51]))
			p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[51:53]) == 1})
			p.commit(ddGdb)
			
			p = Prediction(ExperimentID, PredictionSet, EGAD_ID, parseFloat(line[53:63]))
			p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[63:65]) == 1})
			p.commit(ddGdb)
			
			p = Prediction(ExperimentID, PredictionSet, FoldX_ID, parseFloat(line[65:75]))
			p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[75:77]) == 1})
			p.commit(ddGdb)
			
			p = Prediction(ExperimentID, PredictionSet, Hunter_ID, parseFloat(line[77:87]))
			p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[87:89]) == 1})
			p.commit(ddGdb)
			
			p = Prediction(ExperimentID, PredictionSet, IMutant2_ID, parseFloat(line[89:99]))
			p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[99:101]) == 1})
			p.commit(ddGdb)
			
			p = Prediction(ExperimentID, PredictionSet, Rosetta_ID, parseFloat(line[101:111]))
			p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[111:113]) == 1})
			p.commit(ddGdb)
			
	F.close()

def parseProTherm(ddGdb):
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
			
			pdbID = record["PDB_wild"][0].upper()
			
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
			
			protein = None
			source = None
			if len(record["PROTEIN"]) >= 1:
				protein = join(record["PROTEIN"], " ")
			if len(record["SOURCE"]) >= 1:
				source = join(record["SOURCE"], " ")
			
			if store:
				# Parse mutant
				mutantlist = []
				if record["PDB_mutant"]:
					mutantlist = record["PDB_mutant"]
					for k in range(len(mutantlist)):
						m = mutantlist[k]
						if len(m) > 4:
							if len(m) == 5 and m[4] == ",":
								mutantlist[k] = m[0:4]
							else:
								raise Exception('Error parsing mutant list "%s" in record %d: ' % (mutantlist, ID))
						AllPDBIDs[mutantlist[k]] = True

				
				# NOTE: THESE ARE PATCHES FOR BOTH BAD AND MISMATCHED DATA IN ProTherm 
				if pdbID in ["1CSP", "2LZM"]:
					chainID = "A"
				if ID in [13535]:
					mutantlist.remove("166H")
				
				s = '''
				Number of mutations: 23359
				Number of acceptable mutations: 2661
				Number of unique mutations: 1632
				Number of potentially correctable single errors: 111
				Number of unique PDB IDs: 63
				Chain IDs: -, 1, A, B, I, _
				Number of mutations with standard deviation > 0.50: 111
'''
				s=''' <= 2.0 resolution
				Number of mutations: 23359
				Number of acceptable mutations: 1814
				Number of unique mutations: 1201
				Number of potentially correctable single errors: 111
				Number of unique PDB IDs: 45
				Chain IDs: -, 1, A, B, I, _
				Number of mutations with standard deviation > 0.50: 79
				'''
				Structure = PDBStructure(pdbID, protein = protein, source = source)
				Structure.commit(ddGdb)
				
				ExperimentID = None#Experiment.getDBID() todo
				
										
				# All necessary information is present. Store the record.
				mutationID = "%s-%s-%s%d%s" % (pdbID, chainID, WildTypeAA, ResidueID, MutantAA)
				existingExperiment = experiments.get(mutationID)
				if existingExperiment:
					for mutant in mutantlist:
						existingExperiment.addMutant(mutant)
					if existingExperiment.getChains() != [chainID]:
						raise Exception("Chain mismatch on reading ProTherm database (we are currently only parsing for once chain per experiment): %s, %s " % (existingExperiment.getChains(), chainID))
					Experiment.addExperimentalScore(ID, ddG)
				else:
					Experiment = ExperimentSet(pdbID, "ProTherm-2008-09-08-23581")
					for mutant in mutantlist:
						Experiment.addMutant(m)
					Experiment.addMutation(chainID, ResidueID, WildTypeAA, MutantAA)
					Experiment.addChain(chainID)
					Experiment.addExperimentalScore(ID, ddG)
					experiments[mutationID] = Experiment
				count += 1
				newlist.append(ID)
		else:
			break	

	highvariancecount = 0
	for mutationID, e in experiments.iteritems():
		e.mergeScores(maxStdDeviation = MAX_STANDARD_DEVIATION)
		if not(e.isEligible()):
			highvariancecount += 1

	for mutationID, e in experiments.iteritems():
		if e.isEligible():
			e.commit(ddGdb)
	
	PDBIDs = [id[0:4] for id in experiments.keys()]
	pd = {}
	for x in PDBIDs:
		AllPDBIDs[x] = True
		pd[x] = True
	
	print("")
	print("Number of mutations: %d" % totalcount)
	print("Number of acceptable mutations: %d" % count)
	print("Number of unique mutations: %d" % len(experiments))
	print("Number of potentially correctable single errors: %d" % singleErrors)
	print("Number of unique PDB IDs: %d" % len(pd))
	print("Chain IDs: %s" % join(sorted(chains.keys()), ", "))
	print("Number of mutations with standard deviation > %0.2f: %d" % (MAX_STANDARD_DEVIATION, highvariancecount))

	#print("listOfNoIDs",listOfNoIDs)
	errors = mutantProcessingErrors + ddGProcessingErrors + chainProcessingErrors 
	for e in errors:
		print(e)
		

def parseUniProt(ddGdb):
	uniprot = os.path.join("rawdata", "uniprotmapping.csv")
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

def dumpPDBIDs():
	F = open("PDBIDs.txt", "w")
	for id in AllPDBIDs.keys():
		if len(id) > 4:
			print(id)
	pdbids = sorted(list(set(db.UniqueIDs.keys()).union(set(AllPDBIDs.keys()))))
	F.write(join(pdbids, "\n"))
	F.close()

def parseRawData():
	ddGdb = db.ddGDatabase()
	parseSensDatabase(ddGdb)
	parsePotapov(ddGdb)
	#parseUniProt(ddGdb)
	parseProTherm(ddGdb)
	dumpPDBIDs()
	
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
