#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
from string import join
import math
import pickle
import re
sys.path.insert(0, "common")
import ddgproject
import colortext

ExperimentSet = ddgproject.ExperimentSet
Prediction = ddgproject.Prediction
PDBStructure = ddgproject.PDBStructure
fn = ddgproject.FieldNames()

AllPDBIDs = {}

MAX_RESOLUTION = 2.1
MAX_NUMRES_PROTHERM = 350
MAX_STANDARD_DEVIATION = 1.0

def parseFloat(str):
	if str.strip() == "NaN":
		return None
	else:
		return float(str)

def pdb2toRes(pdbID):
	''' Used for debugging.'''
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
	
def convertSensDataset():
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
				assert(line[4].isdigit())
				resid = int(line[4])
				mut = line[6].strip().split("_")[1]
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
					mutres = line[6].split("_")
					resid = int(mutres[0])
					mut = mutres[1].strip()
					
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
			
def parseSensDataset(ddGdb):
	colortext.message("Parsing Sen's dataset")
	
	mutantsfile = os.path.join("rawdata", "sens-complex-dataset-exp-ddg-reformatted.csv")
	headers = ['Interface', 'Structure', 'Mutant' ,'Chains', 'Mutations', 'ddGexp']
	recordlength = len(headers)
	F = open(mutantsfile)
	linecount = 1
	colortext.printf("|*********************|")			#Progress meter
	for line in F.read().split("\n")[1:]:
		if linecount % 16 == 0:
			colortext.write(".", "lightgreen")
			colortext.flush()
		
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
			
			Experiment = ExperimentSet(pdbID, "SenLiu-ComplexExperimentalDataset", interface = interface)
			if mutant:
				Experiment.addMutant(mutant)
			for mutation in mutations:
				chainID = mutation[0]
				wildtypeAA = mutation[2]
				residueID = mutation[3:-1]
				mutantAA = mutation[-1]
				Experiment.addMutation(chainID, residueID, wildtypeAA, mutantAA, ID = linecount)
			for chainID in chainIDs:
				 Experiment.addChain(chainID)
			Experiment.addExperimentalScore(linecount, ddG, pdbID)
			ExperimentID = Experiment.commit(ddGdb)
			
			linecount += 1
	colortext.printf("")
	F.close()
	
def parsePotapov(ddGdb):
	colortext.message("Parsing Potapov")
	
	PredictionSet = "Potapov-2009"
	
	CCPBSA_ID = ddGdb.execute("SELECT ID FROM Protocol WHERE ID='Potapov:10.1093/protein/gzp030::CC/PBSA';")
	assert(len(CCPBSA_ID) == 1)
	CCPBSA_ID = CCPBSA_ID[0]["ID"]
	
	EGAD_ID = ddGdb.execute("SELECT ID FROM Protocol WHERE ID='Potapov:10.1093/protein/gzp030::EGAD';")
	assert(len(EGAD_ID) == 1)
	EGAD_ID = EGAD_ID[0]["ID"]
	
	FoldX_ID = ddGdb.execute("SELECT ID FROM Protocol WHERE ID='Potapov:10.1093/protein/gzp030::FoldXv3.0';")
	assert(len(FoldX_ID) == 1)
	FoldX_ID = FoldX_ID[0]["ID"]
	
	Hunter_ID = ddGdb.execute("SELECT ID FROM Protocol WHERE ID='Potapov:10.1093/protein/gzp030::Hunter';")
	assert(len(Hunter_ID) == 1)
	Hunter_ID = Hunter_ID[0]["ID"]
	
	IMutant2_ID = ddGdb.execute("SELECT ID FROM Protocol WHERE ID='Potapov:10.1093/protein/gzp030::IMutant2v2.0';")
	assert(len(IMutant2_ID) == 1)
	IMutant2_ID = IMutant2_ID[0]["ID"]
	
	Rosetta_ID = ddGdb.execute("SELECT ID FROM Protocol WHERE ID='Potapov:10.1093/protein/gzp030::Rosettav2.2.0';")
	assert(len(Rosetta_ID) == 1)
	Rosetta_ID = Rosetta_ID[0]["ID"]
	
	results = ddGdb.execute(("SELECT * FROM Experiment WHERE %s" % fn.Source) + "=%s;", parameters = (PredictionSet,))
	if len(results) != 10:
		mutantsfile = os.path.join("rawdata", "mutants.txt")
		F = open(mutantsfile)
		
		colortext.printf("|*********************|")			#Progress meter
		count = 1
		for line in F.read().split("\n"):
			if count % 93 == 0:
				colortext.write(".", "lightgreen")
				colortext.flush()
			count += 1
			
			if line.strip() and line[0] != '#':
				chainID = line[9]
				pdbID = line[5:9].upper()
				AllPDBIDs[pdbID] = True
				
				Structure = PDBStructure(pdbID)
				Structure.commit(ddGdb)
				
				PotapovID = int(line[0:4])
				Experiment = ExperimentSet(pdbID, "Potapov-2009")
				
				Experiment.addMutation(chainID, line[20:26], line[14], line[19], ID = PotapovID) # Include the insertion code
				
				Experiment.addChain(chainID)
				Experiment.addExperimentalScore(PotapovID, parseFloat(line[26:35]), pdbID, numMeasurements = int(line[36:41]))
				ExperimentID = Experiment.commit(ddGdb)
				
				score = parseFloat(line[41:51]) 
				if score:
					score = pickle.dumps({"type" : "simple", "data" : score})
					p = Prediction(ExperimentID, PredictionSet, CCPBSA_ID, score, 'done')
					p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[51:53]) == 1})
					p.commit(ddGdb)
				
				score = parseFloat(line[53:63])
				if score:
					score = pickle.dumps({"type" : "simple", "data" : score})
					p = Prediction(ExperimentID, PredictionSet, EGAD_ID, score, 'done')
					p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[63:65]) == 1})
					p.commit(ddGdb)
					
				score = parseFloat(line[65:75])
				if score:
					score = pickle.dumps({"type" : "simple", "data" : score})
					p = Prediction(ExperimentID, PredictionSet, FoldX_ID, score, 'done')
					p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[75:77]) == 1})
					p.commit(ddGdb)
					
				score = parseFloat(line[77:87])
				if score:
					score = pickle.dumps({"type" : "simple", "data" : score})
					p = Prediction(ExperimentID, PredictionSet, Hunter_ID, score, 'done')
					p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[87:89]) == 1})
					p.commit(ddGdb)
					
				score = parseFloat(line[89:99])
				if score:
					score = pickle.dumps({"type" : "simple", "data" : score})
					p = Prediction(ExperimentID, PredictionSet, IMutant2_ID, score, 'done')
					p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[99:101]) == 1})
					p.commit(ddGdb)
					
				score = parseFloat(line[101:111])
				if score:
					score = pickle.dumps({"type" : "simple", "data" : score})
					p = Prediction(ExperimentID, PredictionSet, Rosetta_ID, score, 'done')
					p.setOptional(Description = {"MutationUsedInEvaluatingTheMethod" : int(line[111:113]) == 1})
					p.commit(ddGdb)
		colortext.printf("")
		F.close()

def parseProTherm(ddGdb):
	'''todo: I realized after the fact that the regexes below do not deal properly with insertion codes in
	   the mutation e.g. "MUTATION        Y 27D D, S 29 D" in record 5438. However, none of these records
	   have ddG values for ProTherm on 2008-09-08 (23581 entries) so we can ignore this issue unless the
	   database is updated.'''

	colortext.message("Parsing ProTherm")
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
	chains = {}
	singleErrors = {} 
	patchthis = {}
	patchPDBs = {}
	
	# These are the records where the mutations include insertion codes
	# It turns out that none of these are eligible for inclusion in the database
	# as they are all missing ddG values.
	# As a precaution, we cast all residue IDs to int which will raise an exception
	# if an insertion code is parsed.
	iCodeRecords = [5438, 5439, 5440, 5441, 8060, 13083, 13084]

	# NOTE: THESE ARE PATCHES FOR MISSING DATA IN ProTherm
	patch = {
		2396 : {'PDB_wild' : None}, # -> 2405. P08505 No related PDB entry. 
		2397 : {'PDB_wild' : None}, #P08505
		2398 : {'PDB_wild' : None}, #P08505
		2400 : {'PDB_wild' : None}, #P08505
		2401 : {'PDB_wild' : None}, #P08505
		2403 : {'PDB_wild' : None}, #P08505
		2404 : {'PDB_wild' : None}, #P08505
		2405 : {'PDB_wild' : None}, #P08505
		4216 : {'PDB_wild' : None}, #P00912 No related PDB entry.
		8588 : {'LENGTH' : 104}, #1ONC
		8589 : {'LENGTH' : 104}, #1ONC
		8590 : {'LENGTH' : 104}, #1ONC
		8591 : {'LENGTH' : 104}, #1ONC
		8592 : {'LENGTH' : 104}, #1ONC
		8593 : {'LENGTH' : 104}, #1ONC
		8594 : {'LENGTH' : 104}, #1ONC
		8595 : {'LENGTH' : 104}, #1ONC
		8596 : {'LENGTH' : 104}, #1ONC
		14229 : {'PDB_wild' : None}, # -> 14233. P08821 No related PDB entry.
		14230 : {'PDB_wild' : None}, #P08821
		14231 : {'PDB_wild' : None}, #P08821
		14232 : {'PDB_wild' : None}, #P08821
		14233 : {'PDB_wild' : None}, #P08821
		14978 : {'LENGTH' : 238}, #1CHK
		14979 : {'LENGTH' : 238}, #1CHK
		14980 : {'LENGTH' : 238}, #1CHK
		14981 : {'LENGTH' : 238}, #1CHK
		14987 : {'LENGTH' : 238}, #1CHK
		14988 : {'LENGTH' : 238}, #1CHK
		14989 : {'LENGTH' : 238}, #1CHK
		14990 : {'LENGTH' : 238}, #1CHK
		14996 : {'LENGTH' : 238}, #1CHK
		14997 : {'LENGTH' : 238}, #1CHK
		14998 : {'LENGTH' : 238}, #1CHK
		14999 : {'LENGTH' : 238}, #1CHK
		16597 : {'LENGTH' : 435}, #1KFW
		16598 : {'LENGTH' : 435}, #1KFW
		16599 : {'LENGTH' : 435}, #1KFW
		16600 : {'LENGTH' : 435}, #1KFW
		19423 : {'MUTATED_CHAIN' : None},# -> 19538. 1OTR A33 - I cannot determine what chain this is
		19424 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19425 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19426 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19427 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19428 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19429 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19430 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19431 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19432 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19433 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19434 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19449 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19450 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19451 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19452 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19453 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19454 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19455 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19456 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19457 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19458 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19459 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19460 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19475 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19476 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19477 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19478 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19479 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19480 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19481 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19482 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19483 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19484 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19485 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19486 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19501 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19502 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19503 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19504 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19505 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19506 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19507 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19508 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19509 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19510 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19511 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19512 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19527 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19528 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19529 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19530 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19531 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19532 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19533 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19534 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19535 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19536 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19537 : {'MUTATED_CHAIN' : None},# 1OTR A33
		19538 : {'MUTATED_CHAIN' : None},# 1OTR A33
		21040 : {'MUTATED_CHAIN' : None},# -> 21332. 1CSP Cannot determine what the mutation is
		21041 : {'MUTATED_CHAIN' : None},# 1CSP
		21097 : {'MUTATED_CHAIN' : None},# 1CSP
		21098 : {'MUTATED_CHAIN' : None},# 1CSP
		21157 : {'MUTATED_CHAIN' : None},# 1CSP
		21158 : {'MUTATED_CHAIN' : None},# 1CSP
		21215 : {'MUTATED_CHAIN' : None},# 1CSP
		21216 : {'MUTATED_CHAIN' : None},# 1CSP
		21273 : {'MUTATED_CHAIN' : None},# 1CSP
		21274 : {'MUTATED_CHAIN' : None},# 1CSP
		21331 : {'MUTATED_CHAIN' : None},# 1CSP
		21332 : {'MUTATED_CHAIN' : None},# 1CSP
	}
	
	# These PDB files have exactly one chain but ProTherm lists the wrong chain e.g. '-' rather than 'A'
	singleChainPDBs = {
		'1A23' : {'MUTATED_CHAIN' : 'A'},
		'1AG2' : {'MUTATED_CHAIN' : 'A'},
		'1AKK' : {'MUTATED_CHAIN' : 'A'},
		'1B5M' : {'MUTATED_CHAIN' : 'A'},
		'1BCX' : {'MUTATED_CHAIN' : 'A'},
		'1BPI' : {'MUTATED_CHAIN' : 'A'},
		'1BTA' : {'MUTATED_CHAIN' : 'A'},
		'1BVC' : {'MUTATED_CHAIN' : 'A'},
		'1CAH' : {'MUTATED_CHAIN' : 'A'},
		'1CSP' : {'MUTATED_CHAIN' : 'A'},
		'1CYO' : {'MUTATED_CHAIN' : 'A'},
		'1FLV' : {'MUTATED_CHAIN' : 'A'},
		'1FTG' : {'MUTATED_CHAIN' : 'A'},
		'1HME' : {'MUTATED_CHAIN' : 'A'},
		'1IOB' : {'MUTATED_CHAIN' : 'A'},
		'1IRO' : {'MUTATED_CHAIN' : 'A'},
		'1L63' : {'MUTATED_CHAIN' : 'A'},
		'1LZ1' : {'MUTATED_CHAIN' : 'A'},
		'1MGR' : {'MUTATED_CHAIN' : 'A'},
		'1ONC' : {'MUTATED_CHAIN' : 'A'},
		'1POH' : {'MUTATED_CHAIN' : 'A'},
		'1RRO' : {'MUTATED_CHAIN' : 'A'},
		'1RTB' : {'MUTATED_CHAIN' : 'A'},
		'1RX4' : {'MUTATED_CHAIN' : 'A'},
		'1SSO' : {'MUTATED_CHAIN' : 'A'},
		'1STN' : {'MUTATED_CHAIN' : 'A'},
		'1SUP' : {'MUTATED_CHAIN' : 'A'},
		'1VQB' : {'MUTATED_CHAIN' : 'A'},
		'1YCC' : {'MUTATED_CHAIN' : 'A'},
		'2ABD' : {'MUTATED_CHAIN' : 'A'},
		'2ACY' : {'MUTATED_CHAIN' : 'A'},
		'2AKY' : {'MUTATED_CHAIN' : 'A'},
		'2HPR' : {'MUTATED_CHAIN' : 'A'},
		'2LZM' : {'MUTATED_CHAIN' : 'A'},
		'2RN2' : {'MUTATED_CHAIN' : 'A'},
		'3SSI' : {'MUTATED_CHAIN' : 'A'},
		'451C' : {'MUTATED_CHAIN' : 'A'},
		'4LYZ' : {'MUTATED_CHAIN' : 'A'},
	}
	
	# In these cases, the data in ProTherm is incorrect according to the publication
	overridden = {
		# These cases have ambiguous entries
		
		
		3469  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1G6N'}, # Two identical chains A and B but '-' specified in ProTherm
		3470  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1G6N'}, 
		14153 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1N0J'}, # Two identical chains A and B but '-' specified in ProTherm
		2418  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'}, # Two identical chains A and B but '-' specified in ProTherm
		5979  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
		5980  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
		5981  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
		5982  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
		5983  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
		5984  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
		5985  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
		5986  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
		5987  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},	
		3629  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'}, # Two identical chains A and B but '-' specified in ProTherm
		3630  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3631  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3632  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3633  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3634  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3635  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3636  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3637  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3638  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3639  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3640  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3641  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3642  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3643  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3644  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
		3604  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'}, # Three identical chains A, B, and C but '-' specified in ProTherm
		3605  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		3606  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		3607  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		3608  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		3609  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		3610  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		3611  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13412 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13413 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13414 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13415 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13416 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13417 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13418 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13419 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13420 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13985 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13986 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		13421 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		14253 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		14254 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		14255 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
		8302  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'}, # Four identical chains A, B, C, and D but '-' specified in ProTherm
		8303  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
		8304  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
		8305  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
		8306  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
		14474 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
		14475 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
		14476 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
		10057 : {'MUTATED_CHAIN' : 'B', 'PDB' : '1RN1'}, # Three identical chains A, B, C. These cases are odd ones - the PDB file contains identical chains but residue 45 is missing in chain A
		10058 : {'MUTATED_CHAIN' : 'B', 'PDB' : '1RN1'},
	}

	# In these cases, the protein is elongated at the 67th position. This is more than a mutation so I ignore it. 	
	skipTheseCases = [12156, 12159, 12161, 12188, 12191, 12193, 14468]
	
	# In these cases, the structural information needed for the mutations (residues A57, A62) is missing in the PDB file
	# Some of the information is present in chain B (A57's corresponding residue) but I decided against remapping the residues. 
	skipTheseCases.extend([13451, 13452])
	
	# In this case, I couldn't map the paper's structure well enough onto 1YCC. The mutations mentioned e.g. A57 N->I do not correspond with the PDB structure (attributed to a later paper).
	skipTheseCases.append(11817)
		
	mutationregex = re.compile("^(\w\d+\w,)+(\w\d+\w)[,]{0,1}$")

	CysteineMutationCases = [13663, 13664, 13677, 13678]	
	multimapCases1 = [22383, 22384]
	mmapCases1regex = re.compile("PDB:(.+[,]{0,1})+\)")
	multimapCases2 = range(17681, 17687 + 1) + range(17692, 17698 + 1)
	mmapCases2regex = re.compile("^.*PDB:(.*);PIR.*$")			
	multimapCases3 = range(14215, 14223 + 1) + [16991, 17678, 17679, 17680, 17689, 17690, 17691]
	mmapCases3regex = re.compile("PDB:(\w.*?\w)[;)]")
	
	# These fields of ProTherm records cannot be empty for our purposes
	requiredFields = ["NO.", "PDB_wild", "LENGTH", "ddG", "MUTATION", "MUTATED_CHAIN"]
	
	failthis = False
	newlist = []
	protherm = os.path.join("rawdata", "ProTherm.dat")
	F = open(protherm)
	colortext.printf("|*********************|")
	while(True):
		# Read a record
		record = {}
		line = F.readline()
		while line and not(line.startswith("//")):
			if line[0] != "*" and line[0] != " ":
				line = line.split()
				if not singleErrors.get(line[0]):
					singleErrors[line[0]] = 0
				if len(line) > 1:
					record[line[0]] = line[1:]
				else:
				 	record[line[0]] = None
			line = F.readline()
					
		# Parse the results
		if record:
			MutantAA = []
			WildTypeAA = []
			ResidueID = []
			chainID = None
			
			totalcount += 1	
			
			# Find out whether we have enough information
			store = True
			ID = int(record["NO."][0])
			
			if overridden.get(ID):
				#print("Overriding case %d" % ID)
				if record.get('PDB_wild'):
					if record["PDB_wild"][0] != overridden[ID]["PDB"]:
						raise colortext.Exception("Error in overridden table: Record %d. Read '%s' for PDB_wild, expected '%s'." % (ID, record["PDB_wild"], overridden[ID]["PDB"]))
				for k,v in overridden[ID].iteritems():
					if k != "PDB":
						record[k] = v.split()
			if record["PDB_wild"]:
				pdbID = record["PDB_wild"][0].upper()
				if singleChainPDBs.get(pdbID):
					for k,v in singleChainPDBs[pdbID].iteritems():
						record[k] = v.split()
				
			#Progress meter
			if ID % 1000 == 0:
				colortext.write(".", "green")
				colortext.flush()
			
			missingFields = []
			if ID in skipTheseCases:
				continue
			for field in requiredFields:
				if not record[field]:
					store = False
					missingFields.append(field)
			if not store:
				if len(missingFields) == 1:
					if record["MUTATION"][0] and record["MUTATION"][0] != "wild" and missingFields[0] != "ddG":
						if not record["MUTATED_CHAIN"]:
							if not patch.get(ID):
								colortext.error("Error processing chain: ID %d, no chain" %  (ID))
								singleErrors["MUTATED_CHAIN"] += 1
								patchthis[ID] = "MUTATED_CHAIN %s-%s" % (record.get("PDB_wild"), record.get("MUTATION")) 
							elif patch[ID]["MUTATED_CHAIN"]:
								record["MUTATED_CHAIN"] = [patch[ID]["MUTATED_CHAIN"]]
								store = True
						elif not record["LENGTH"]:
							if not patch.get(ID):
								colortext.error("Error processing length: ID %d, no length" %  (ID))
								singleErrors["LENGTH"] += 1
								patchthis[ID] = "LENGTH %s" % record.get("PDB_wild") 
							elif patch[ID]["LENGTH"]:
								record["LENGTH"] = [patch[ID]["LENGTH"]]
								store = True
						elif not record["PDB_wild"]:
							if not patch.get(ID):
								colortext.error("Error processing PDB ID: ID %d, no PDB ID" %  (ID))
								singleErrors["PDB_wild"] += 1
								patchthis[ID] = "PDB_wild: %s" % record.get("SWISSPROT_ID") 
							elif patch[ID]["PDB_wild"]:
								record["PDB_wild"] = [patch[ID]["PDB_wild"]]
								store = True
						else:
							colortext.error("Error processing structure: ID %d, no %s " % (ID, missingFields[0]))
							singleErrors[missingFields[0]] += 1
				if not store:
					continue
		
			# Only allow proteins with <= MAX_NUMRES_PROTHERM residues
			if int(record["LENGTH"][0]) > MAX_NUMRES_PROTHERM:
				continue
		
			pdbID = record["PDB_wild"][0].upper()
			
			# Parse chain
			if len(record["MUTATED_CHAIN"]) == 1:
				chainID = record["MUTATED_CHAIN"][0]
				chains[chainID] = True
			else:
				colortext.error("Error processing chain: ID %d, %s" %  (ID, record["MUTATED_CHAIN"]))
				store = False
			
			# Parse mutation
			mutationline = record["MUTATION"]
			if ID in CysteineMutationCases: # Hack for input which includes information on the bonds generated from mutation to cysteine
				ResidueID = [int(mutationline[1])] 
				WildTypeAA = [mutationline[0]]
				MutantAA = [mutationline[2]]
			elif ID in multimapCases1:
				cline = join(mutationline, "")
				mtchslst = mmapCases1regex.findall(cline)
				if mtchslst:
					assert(len(mtchslst) == 1)
					mtchslst = mtchslst[0].split(',')
					for mtch in mtchslst:
						assert(mtch)
						ResidueID.append(int(mtch[1:-1]))
						WildTypeAA.append(mtch[0])
						MutantAA.append(mtch[-1])
				else:
					raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
			elif ID in multimapCases2:
				cline = join(mutationline, "")
				mtchslst = mmapCases2regex.match(cline)
				if mtchslst:
					mtchslst = mtchslst.group(1).split(",")
					for mtch in mtchslst:
						assert(mtch)
						ResidueID.append(int(mtch[1:-1]))
						WildTypeAA.append(mtch[0])
						MutantAA.append(mtch[-1])
				else:
					raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
			elif ID in multimapCases3:
				cline = join(mutationline, "")
				mtchslst = mmapCases3regex.findall(cline)
				if mtchslst:
					for mtch in mtchslst:
						assert(mtch)
						ResidueID.append(int(mtch[1:-1]))
						WildTypeAA.append(mtch[0])
						MutantAA.append(mtch[-1])
				else:
					raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
			elif len(mutationline) == 3:
				try:
					ResidueID = [int(mutationline[1])]
					WildTypeAA = [mutationline[0]]
					MutantAA = [mutationline[2]]
				except:
					cline = join(mutationline, "")
					raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
			elif len(mutationline) == 1:
				mline = mutationline[0]
				if mline == "wild" or mline == "wild*" or mline == "wild**":
					store = False
				else:
					cline = join(mutationline, "")
					raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
			elif len(mutationline) % 3 == 0 or len(mutationline) == 5: #2nd case is a hack for 11146
				cline = join(mutationline, "").strip()
				m = mutationregex.match(cline)
				if m:
					mutations = [m.group(i)[:-1] for i in range(1, m.lastindex)] + [m.group(m.lastindex)]
					for mut in mutations:
						ResidueID.append(int(mut[1:-1]))
						WildTypeAA.append(mut[0])
						MutantAA.append(mut[-1])
				else:
					raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
				 
			else:
				cline = join(mutationline, "")
				raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
				
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
						colortext.error("Error processing ddG: ID %d, %s" % (ID, record["ddG"]))
						store = False
				else:
					idx = ddGline.find("kcal/mol")
					if idx != -1:
						try:
							ddG = 4.184 * float(ddGline[0:idx])
						except:
							colortext.error("Error processing ddG: ID %d, %s" % (ID, record["ddG"]))
							store = False
					else:
						colortext.error("Error processing ddG: ID %d, %s" % (ID, record["ddG"]))
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
				if ID in [13535]:
					mutantlist.remove("166H")
				
				Structure = PDBStructure(pdbID, protein = protein, source = source)
				Structure.commit(ddGdb)
				
				# All necessary information is present. Store the record.
				assert(ResidueID and len(ResidueID) == len(WildTypeAA) and len(WildTypeAA) == len(MutantAA))
				mutuationIDsuffix = join(["%s%s%s" % (WildTypeAA[i], ResidueID[i], MutantAA[i]) for i in range(len(WildTypeAA))], "-")
				mutationID = "%s-%s-%s" % (pdbID, chainID, mutuationIDsuffix)
				Experiment = None
				existingExperiment = experiments.get(mutationID)
				if existingExperiment:
					for mutant in mutantlist:
						existingExperiment.addMutant(mutant)
					if existingExperiment.getChains() != [chainID]:
						raise Exception("Chain mismatch on reading ProTherm database (we are currently only parsing for once chain per experiment): %s, %s " % (existingExperiment.getChains(), chainID))
					try:
						existingExperiment.addExperimentalScore(ID, ddG, pdbID)
					except Exception, e:
						colortext.error((ID, mutationID))
						colortext.error(str(e))
						raise Exception()
				else:
					Experiment = ExperimentSet(pdbID, "ProTherm-2008-09-08-23581")
					for mutant in mutantlist:
						Experiment.addMutant(m)
					
					test = []
					for midx in range(len(ResidueID)):
						if ResidueID[midx] not in test:
							test.append(ResidueID[midx])
						else:
							failthis = True
							colortext.error("Error in record %d" % ID)
					for midx in range(len(ResidueID)):
						Experiment.addMutation(chainID, ResidueID[midx], WildTypeAA[midx], MutantAA[midx], ID = ID)
					Experiment.addChain(chainID)
					try:
						Experiment.addExperimentalScore(ID, ddG, pdbID)
					except:
						colortext.error(str(e))
						raise Exception()
					experiments[mutationID] = Experiment
				count += 1
				newlist.append(ID)
		else:
			break	
	colortext.printf("\t[Parsing complete]", 'green')
	
	if failthis:
		sys.exit(0)
	highvariancecount = 0
	for mutationID, e in experiments.iteritems():
		e.mergeScores(maxStdDeviation = MAX_STANDARD_DEVIATION)
		if not(e.isEligible()):
			highvariancecount += 1
	
	PROGRESS_BAR_LENGTH = 23
	progressCounter = 1
	progressStep = (len(experiments) / PROGRESS_BAR_LENGTH) or 2
	for mutationID, e in experiments.iteritems():
		#Progress meter
		if progressCounter % progressStep == 0:
			colortext.write(".", "lightgreen")
			colortext.flush()
		ExperimentID = e.commit(ddGdb)
		progressCounter += 1
	colortext.printf("\t[Database storage complete]", 'lightgreen')
		
	PDBIDs = [id[0:4] for id in experiments.keys()]
	pd = {}
	for x in PDBIDs:
		AllPDBIDs[x] = True
		pd[x] = True
	
	print("")
	print("Number of mutations: %d" % totalcount)
	print("Number of acceptable mutations: %d" % count)
	print("Number of unique mutations: %d" % len(experiments))
	se = 0
	for k, count in singleErrors.iteritems():
		se += count
	print("Number of potentially correctable single errors: %d" % se)
	for k, count in singleErrors.iteritems():
		if k in requiredFields and count != 0:
			print("\t%s: %d" % (k, count))
	print("Number of unique PDB IDs: %d" % len(pd))
	print("Chain IDs: %s" % join(sorted(chains.keys()), ", "))
	print("Number of mutations with standard deviation > %0.2f: %d" % (MAX_STANDARD_DEVIATION, highvariancecount))
	
	if patchthis:
		print("patch = {")
		for id, type in sorted(patchthis.iteritems()):
			print("%d : {'%s' : ''}," % (id, type))
		print("}")
	
	#print("listOfNoIDs",listOfNoIDs)
	errors = mutantProcessingErrors + ddGProcessingErrors + chainProcessingErrors 
	for e in errors:
		print(e)
		
def dumpPDBIDs():
	F = open("PDBIDs.txt", "w")
	for id in AllPDBIDs.keys():
		if len(id) > 4:
			print(id)
	pdbids = sorted(list(set(ddgproject.UniqueIDs.keys()).union(set(AllPDBIDs.keys()))))
	F.write(join(pdbids, "\n"))
	F.close()

def parseExtraProThermData(ddGdb):
	colortext.message("Parsing ProTherm")
	protherm = os.path.join("rawdata", "ProTherm.dat")
	F = open(protherm)
	colortext.printf("|*********************|")
	
	DBIDs = {}
	for r in ddGdb.execute('SELECT SourceID, ExperimentScore.ID AS ID FROM ExperimentScore INNER JOIN Experiment ON ExperimentScore.ExperimentID = Experiment.ID WHERE Experiment.Source="ProTherm-2008-09-08-23581"', cursorClass = ddgproject.StdCursor):
		DBIDs[int(r[0])] = int(r[1])
	count=  0
	
	exp2DBfield = {
		"T"				: fn.ExpConTemperature,
		"pH"			: fn.ExpConpH,
		"BUFFER_NAME"	: fn.ExpConBuffer,
		"BUFFER_CONC"	: fn.ExpConBufferConcentration,
		"ION_NAME_1"	: fn.ExpConIon,
		"ION_CONC_1"	: fn.ExpConIonConcentration,
		"PROTEIN_CONC"	: fn.ExpConProteinConcentration,
		"MEASURE"		: fn.ExpConMeasure,
		"METHOD"		: fn.ExpConMethodOfDenaturation,
	}
	expfields = exp2DBfield.keys()
	numericexpfields = ["T", "pH"]  
	missingExpData = {}
	missingReferences = 0
	noPMIDs = {}
	PMIDlist = {}
	maxDBfieldlengths = {}
	
	expCondSQL = "UPDATE ExperimentScore SET"
	for f in sorted(exp2DBfield.values()):
		expCondSQL += ((" %s=" % f) + "%s,")
	expCondSQL = expCondSQL[:-1]
	expCondSQL += " WHERE ID=%s;"
	
	missingRefMap = {
		"BIOCHEMISTRY 34, 7094-7102 (1995)" 		: ("PMID", 7766619),
		"BIOCHEMISTRY 34, 7103-7110 (1995)" 		: ("PMID", 7766620), # This is probably the correct record. Pubmed has a different end page (two extra pages)
		"BIOCHEMISTRY 37, 2477-2487 (1998)" 		: ("PMID", 9485396),
		"BIOCHIM BIOPHYS ACTA 1429, 365-376 (1999)" : ("PMID", 9989221),
		"PROTEIN SCI 6, 2196-2202 (1997)" 			: ("PMID", 9336842),
	} 
	
	publicationSources = {}
	for r in ddGdb.execute('SELECT ID FROM %s' % fn.Source, cursorClass = ddgproject.StdCursor):
		publicationSources[r[0]] = True
	
	count = 0
	while(True):
		# Read a record
		record = {}
		line = F.readline()
		while line and not(line.startswith("//")):
			if line[0] != "*" and line[0] != " ":
				line = line.split()
				#if not singleErrors.get(line[0]):
				#	singleErrors[line[0]] = 0
				if len(line) > 1:
					record[line[0]] = line[1:]
				else:
				 	record[line[0]] = None
			line = F.readline()
		
		
		kelvin_regex = re.compile("^(\d*[.]?\d*) K$")
		celsius_regex = re.compile("^(\d*[.]?\d*) C$")
							
		# Parse the results
		if record:
			# Find out whether we have enough information
			store = True
			ID = int(record["NO."][0])
			if DBIDs.get(ID):
				ExperimentScoreID = DBIDs[ID]
				count +=1
				
				ExperimentalScoreUpdate = {}
				
				# Add experimental condition 
				for h in expfields:
					if not record[h]:
						missingExpData[h] = missingExpData.get(h, 0) + 1
					fielddata = record[h]
					if fielddata:
						fielddata = join(fielddata, " ")
					if record[h]:
						maxDBfieldlengths[h] = max(maxDBfieldlengths.get(h, 0), len(fielddata))
					if fielddata == "Unknown":
						fielddata = None
					if fielddata and h in numericexpfields:
						try:
							fielddata = float(fielddata)
						except:
							if h == "T":
								K = kelvin_regex.match(fielddata)
								C = celsius_regex.match(fielddata)
								if K:
									fielddata = float(K.group(1)) - 273.15 
								elif C:
									fielddata = float(C.group(1))
								else:
									fielddata = None
						if fielddata == None:
							colortext.error("Failed to convert field %s's number %s for record %d." % (h, fielddata, ID))
							
					ExperimentalScoreUpdate[exp2DBfield[h]] = fielddata or None
				
				updateVals = tuple([ExperimentalScoreUpdate[f] for f in sorted(exp2DBfield.values())] + [ExperimentScoreID]) 
				ddGdb.execute(expCondSQL, parameters = updateVals)
				
				# Add reference
				if not record["REFERENCE"]:
					missingReferences += 1
				referenceData = record["REFERENCE"] or ""
				referenceData = join(referenceData, " ").strip()
				idx = referenceData.find("PMID:")
				referenceID = None
				if idx != -1:
					refPMID = referenceData[idx+5:].strip()
					if refPMID:
						if refPMID.isdigit():
							referenceID = int(refPMID)
						else:
							colortext.error("Check reference for record %(ID)d: %(referenceData)s. It is not numeric." % vars())				
				if not referenceID:
					refdetails = None
					if idx == -1:
						refdetails = referenceData
					else:
						refdetails = referenceData[:idx]
					if missingRefMap.get(refdetails) and missingRefMap[refdetails][0] == "PMID":
						referenceID = missingRefMap[refdetails][1]
					else:
						missingReferences += 1
						authors = record["AUTHOR"]
						colortext.warning("No PMID reference for record %(ID)d: '%(referenceData)s', %(authors)s." % vars())
						noPMIDs[referenceData.strip()] = authors
				
				if referenceID:
					assert(str(referenceID).isdigit())
					dbReferencePK = "PMID:%s" % referenceID
					if not publicationSources.get(dbReferencePK):
						print("Adding source '%s'." % dbReferencePK)
						ddGdb.insertDict(fn.Source, {fn.ID : dbReferencePK})
						publicationSources[dbReferencePK] = True
					results = ddGdb.execute("SELECT * FROM SourceLocation WHERE SourceID=%s AND Type='PMID'", parameters = (dbReferencePK,))
					if results:
						assert(len(results) == 1)
						result = results[0]
						if result["ID"] != str(referenceID):
							colortext.error("Source %s in the database has the PMID '%s' whereas we retrieved a value of '%s' from ProTherm." % (dbReferencePK, result["ID"], referenceID))
					else:
						ddGdb.insertDict(fn.SourceLocation, {
							fn.SourceID	: dbReferencePK,
							fn.ID		: referenceID,
							fn.Type		: "PMID",
						})
					
					ddGdb.execute("UPDATE ExperimentScore SET Publication=%s WHERE ID=%s", parameters = (dbReferencePK, int(ExperimentScoreID)))
					

				
			#Progress meter
			if ID % 1000 == 0:
				colortext.write(".", "green")
				colortext.flush()
			
			
		else:
			break
	
	F.close()
	colortext.printf("\n[Parsing complete]\n\n", 'green')
	colortext.printf("Summary", 'green')
	colortext.printf("*******\n", 'green')
	colortext.message("Records counted: %d" % count)
	colortext.message("Number of unique references (PMIDs): %d\n" % len(PMIDlist))
	colortext.message("Field lengths: %s" % maxDBfieldlengths)
	if missingExpData:
		colortext.warning("Missing experimental data")
		for h, v in missingExpData.iteritems():
			colortext.warning("\t%(h)s: for %(v)s records" % vars())
	if missingReferences:
		colortext.warning("Missing references for %d records." % missingReferences)
		for ref, authors in sorted(noPMIDs.iteritems()):
			colortext.warning("\t%s: %s" % (ref, authors))
	
def parseRawData(ddGdb):
	print("")
	#parsePotapov(ddGdb)
	#parseSensDataset(ddGdb)
	#parseProTherm(ddGdb)
	parseExtraProThermData(ddGdb)
	return
	if ddGdb.chainWarnings:
		colortext.warning("UNAMBIGUOUS BADLY-SPECIFIED MUTATIONS")
		for pdbID, recordsAndChain in sorted(ddGdb.chainWarnings.iteritems()):
			colortext.warning("%s: Chains %s" % (pdbID, join(ddgproject.PDBChains.get(pdbID, ["Missing data"]), ", ")))
			for rc in recordsAndChain:
				break
				records = rc[0]
				chain = rc[1]
				colortext.warning("\tRecords: %s, chain %s" % (join(map(str, records), ", "), chain))
	if ddGdb.chainErrors:
		colortext.error("AMBIGUOUS BADLY-SPECIFIED MUTATIONS")
		for pdbID, recordsAndChain in sorted(ddGdb.chainErrors.iteritems()):
			colortext.error("%s: Chains %s" % (pdbID, join(ddgproject.PDBChains.get(pdbID, ["Missing data"]), ", ")))
			for rc in recordsAndChain:
				records = rc[0]
				chain = rc[1]
				colortext.error("\tRecords: %s, chain %s" % (join(map(str, records), ", "), chain))
	dumpPDBIDs()

def main():	
	parseRawData(ddgproject.ddGDatabase())
	
if __name__ == "__main__":
	#print("Preventing accidental runs. Exiting.")
	#sys.exit(0)
	main()
