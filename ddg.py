#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
from string import join
import math
import re
import time
import subprocess
from webparser import MAX_RESOLUTION, MAX_NUMRES_PROTHERM, MAX_STANDARD_DEVIATION
sys.path.insert(0, "common")
import ddgproject
from pdb import PDB, ResidueID2String, checkPDBAgainstMutations, aa1
from Bio.PDB import *
import colortext
import traceback
import pickle


#Rosetta3.1 release. mini revision 30964. minirosetta_database revision 30967.
#Rosetta3.2 release. mini revision 39284
#Rosetta3.3 release: mini revision 42942

# MySQL 5.0 has an inefficient stored procedure implementation for non-persistent connections
# so we will define most SQL queries here.
SQLQueries = {
	"SmallToLarge"				:	'SELECT ExperimentID, Structure, Source, Chain, ResidueID, WildtypeAA, MutantAA, wildtype.Size, mutant.Size FROM Experiment INNER JOIN ExperimentMutation ON ID=ExperimentID INNER JOIN AminoAcid as wildtype ON WildtypeAA = wildtype.Code INNER JOIN AminoAcid as mutant ON MutantAA = mutant.Code WHERE wildtype.Size="small" AND mutant.Size="large";',
	"SmallToLargeD"				:	'SELECT DISTINCT ExperimentID FROM Experiment INNER JOIN ExperimentMutation ON ID=ExperimentID INNER JOIN AminoAcid as wildtype ON WildtypeAA = wildtype.Code INNER JOIN AminoAcid as mutant ON MutantAA = mutant.Code WHERE wildtype.Size="small" AND mutant.Size="large";',
	"LargeToSmall"				:	'SELECT ExperimentID, Structure, Source, Chain, ResidueID, WildtypeAA, MutantAA, wildtype.Size, mutant.Size FROM Experiment INNER JOIN ExperimentMutation ON ID=ExperimentID INNER JOIN AminoAcid as wildtype ON WildtypeAA = wildtype.Code INNER JOIN AminoAcid as mutant ON MutantAA = mutant.Code WHERE wildtype.Size="large" AND mutant.Size="small";',
	"LargeToSmallD"				:	'SELECT DISTINCT ExperimentID FROM Experiment INNER JOIN ExperimentMutation ON ID=ExperimentID INNER JOIN AminoAcid as wildtype ON WildtypeAA = wildtype.Code INNER JOIN AminoAcid as mutant ON MutantAA = mutant.Code WHERE wildtype.Size="large" AND mutant.Size="small";',
	"SmallToSmall"				:	'SELECT ExperimentID, Structure, Source, Chain, ResidueID, WildtypeAA, MutantAA, wildtype.Size, mutant.Size FROM Experiment INNER JOIN ExperimentMutation ON ID=ExperimentID INNER JOIN AminoAcid as wildtype ON WildtypeAA = wildtype.Code INNER JOIN AminoAcid as mutant ON MutantAA = mutant.Code WHERE wildtype.Size="small" AND mutant.Size="small";',
	"LargeToLarge"				:	'SELECT ExperimentID, Structure, Source, Chain, ResidueID, WildtypeAA, MutantAA, wildtype.Size, mutant.Size FROM Experiment INNER JOIN ExperimentMutation ON ID=ExperimentID INNER JOIN AminoAcid as wildtype ON WildtypeAA = wildtype.Code INNER JOIN AminoAcid as mutant ON MutantAA = mutant.Code WHERE wildtype.Size="large" AND mutant.Size="large";',
	"ProThermWithinResolution"	:	'SELECT Experiment.*, Structure.Techniques, Structure.Resolution FROM Experiment INNER JOIN Structure ON Experiment.Structure = PDB_ID WHERE Experiment.Source = "ProTherm-2008-09-08-23581" AND Structure.Resolution <= %f;' % MAX_RESOLUTION,
	"AllProTherm"				:	'SELECT Experiment.*, Structure.Techniques, Structure.Resolution FROM Experiment INNER JOIN Structure ON Experiment.Structure = PDB_ID WHERE Experiment.Source = "ProTherm-2008-09-08-23581";',
	}

StoredProcedures = ["GetScores", "GetChains", "GetInterfaces", "GetMutants", "GetMutations"]
dbfields = ddgproject.FieldNames()

# todo - for the daemon
class JobTestRunner(object):
	'''This class is just here to test the running of jobs at an initial stage of development.
		Do not rely on it being maintained.''' 
	
	def __init__(self):
		self.ddGdb = ddgproject.ddGDatabase()
	
	def __del__(self):
		self.ddGdb.close()
	
	def run(self):
		results = self.ddGdb.execute(("SELECT * FROM %(Prediction)s WHERE Status=" % dbfields) + "%s", parameters = (dbfields.queued,))
		for result_dict in results:
			runPrediction(result_dict)
			
	def runPrediction(self, predictionRecord):
		'''This function takes in a database record with the prediction details and starts
			running the prediction.'''
		ddGdb = self.ddGdb
		try:
			experimentID = predictionRecord[dbfields.ExperimentID]
			parameters = (experimentID,)
			
			chains = [result[0] for result in ddGdb.callproc("GetChains", parameters = parameters, cursorClass = ddgproject.StdCursor)]
			mutants = [result[0] for result in ddGdb.callproc("GetMutants", parameters = parameters, cursorClass = ddgproject.StdCursor)]
			interfaces = [result[0] for result in ddGdb.callproc("GetInterfaces", parameters = parameters, cursorClass = ddgproject.StdCursor)]
			mutations = ddGdb.callproc("GetMutations", parameters = parameters)
		except:
			pass

class JobInserter(object):
	'''This class is responsible for inserting prediction jobs to the database.''' 
	
	def __init__(self):
		self.ddGdb = ddgproject.ddGDatabase()
	
	def __del__(self):
		self.ddGdb.close()
	
	def createResfile(self, pdb, mutations):
		'''The mutations here are in the original PDB numbering. pdb is assumed to use Rosetta numbering.
			We use the pdb mapping from PDB numbering to Rosetta numbering to generate the resfile.
		'''
		resfile = []
		for mutation in mutations:
			chain = mutation[0]
			resid = mutation[1]
			wt = mutation[2]
			mt = mutation[3]
			
			# Check that the expected wildtype exists in the PDB 
			readwt = pdb.getAminoAcid(pdb.getAtomLine(chain, resid))
			assert(wt == aa1[readwt])
			resid = resid.strip()
			resfile.append("%(resid)s %(chain)s PIKAA %(mt)s" % vars())
		if resfile:
			resfile = ["NATAA", "start"] + resfile
			return join(resfile, "\n")
		else:
			raise Exception("An error occurred creating a resfile for the ddG job.")
	
	def add(self, experimentID, PredictionSet, ProtocolID, KeepHETATMLines, Description = {}, InputFiles = {}):
		'''This function inserts a prediction into the database.
			The parameters define:
				the experiment we are running the prediction for;
				the name of the set of predictions for later grouping;
				the short description of the Command to be used for prediction;
				whether HETATM lines are to be kept or not.
			We strip the PDB based on the chains used for the experiment and KeepHETATMLines.
			We then add the prediction record, including the stripped PDB and the inverse mapping
			from Rosetta residue numbering to PDB residue numbering.''' 
			
		parameters = (experimentID,)
		
		try:
			sql = ("SELECT %(PDB_ID)s, %(Content)s FROM %(Experiment)s INNER JOIN %(Structure)s WHERE %(Experiment)s.%(Structure)s=%(PDB_ID)s AND %(Experiment)s.%(ID)s=" % dbfields) + "%s"
			results = self.ddGdb.execute(sql, parameters = parameters)
			if len(results) != 1:
				raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
			
			# Get the related PDB ID and file
			result = results[0]
			pdbID = result[dbfields.PDB_ID]
			contents = result[dbfields.Content]
			
			pdb = PDB(contents.split("\n"))
			
			# Check that the mutated positions exist and that the wild-type matches the PDB
			mutations = [result for result in self.ddGdb.callproc("GetMutations", parameters = parameters, cursorClass = ddgproject.StdCursor)]
			checkPDBAgainstMutations(pdbID, pdb, mutations)
			
			# Strip the PDB to the list of chains. This also renumbers residues in the PDB for Rosetta.
			chains = [result[0] for result in self.ddGdb.callproc("GetChains", parameters = parameters, cursorClass = ddgproject.StdCursor)]
			pdb.stripForDDG(chains, KeepHETATMLines)
			
			# - Post stripping checks -
			# Get the 'Chain ResidueID' PDB-formatted identifier for each mutation mapped to Rosetta numbering
			# then check again that the mutated positions exist and that the wild-type matches the PDB
			remappedMutations = pdb.remapMutations(mutations, pdbID)
			remappedMutations = [[m[0], ResidueID2String(m[1]), m[2], m[3]] for m in remappedMutations]
			
			resfile = self.createResfile(pdb, remappedMutations)
			
			# Check to make sure that we haven't stripped all the ATOM lines
			if not pdb.GetAllATOMLines():
				raise colortext.Exception("No ATOM lines remain in the stripped PDB file of %s." % pdbID)				
			
			# Check to make sure that CSE and MSE are not present in the PDB
			badresidues = pdb.CheckForPresenceOf(["CSE", "MSE"])
			if badresidues:
				raise colortext.Exception("Found residues [%s] in the stripped PDB file of %s. These should be changed to run this job under Rosetta." % (join(badresidues, ", "), pdbID))
			
			# Turn the lines array back into a valid PDB file				
			strippedPDB = join(pdb.lines, "\n")
		except Exception, e:
			colortext.error("\nError: '%s'.\n" % (str(e)))
			colortext.error(traceback.format_exc())
			raise colortext.Exception("An exception occurred retrieving the experimental data for Experiment ID #%s." % experimentID)
		
		InputFiles["RESFILE"] = resfile
		
		ExtraParameters = {}
		InputFiles = pickle.dumps(InputFiles)
		Description = pickle.dumps(Description)
		ExtraParameters = pickle.dumps(ExtraParameters)
		
		params = {
			dbfields.ExperimentID : experimentID,
			dbfields.PredictionSet : PredictionSet,
			dbfields.ProtocolID : ProtocolID,
			dbfields.KeptHETATMLines : KeepHETATMLines,
			dbfields.StrippedPDB : strippedPDB,
			dbfields.ResidueMapping : pickle.dumps(pdb.get_ddGInverseResmap()),
			dbfields.InputFiles : InputFiles,
			dbfields.Description : Description,
			dbfields.Status : dbfields.queued,
			dbfields.ExtraParameters : ExtraParameters,
			#"mutations" : mutations,
			#"remappedmutations" : remappedMutations,
			#"pdbID" : pdbID,
			#dbfields.Command : {
			#	"Preminimization" : [
			#		'%(EXECUTABLE)s',
			#		'-in:file:l', '%(INPUT_PDB_LIST)s',
			#		'-in:file:fullatom',
			#		'-ignore_unrecognized_res',
			#		'-fa_max_dis', '9.0',
			#		'-database', '%(DATABASE)s',
			#		'-ddg::harmonic_ca_tether', '0.5',
			#		'-score:weights', 'standard',
			#		'-ddg::constraint_weight','1.0',
			#		'-ddg::out_pdb_prefix', 'min_cst_0.5',
			#		'-ddg::sc_min_only', 'false',
			#		'-score:patch', 'score12'
			#	]}
		}
		self.ddGdb.insertDict('Prediction', params)
							

def addAllEligibleProTherm( PredictionSet, CommandName, KeepHETATMLines):
	inserter = JobInserter()
	colortext.printf("\nAdding ProTherm mutations to %s prediction set." % PredictionSet, "lightgreen")
	ddGdb = ddgproject.ddGDatabase()
	WithinResolution = ddGdb.execute(SQLQueries["ProThermWithinResolution"])
	experimentIDs = []
	moreThanOneChain = 0
	for exp in WithinResolution:
		experimentID = exp["ID"] 
		if len(ddGdb.callproc("GetMutations", parameters = (experimentID,))) == 1:
			if exp["Techniques"].find('X-RAY DIFFRACTION') != -1:
				if ddGdb.getStandardDeviation(experimentID) <= MAX_STANDARD_DEVIATION:
						if len(ddGdb.callproc("GetChains", parameters = (experimentID,))) == 1:
							experimentIDs.append(experimentID)
						else:
							moreThanOneChain += 1
	colortext.message("\nThe number of unique ProTherm experiments with:\n\t- one mutation;\n\t- structures solved by X-ray diffraction and with <= %d residues;\n\t- a maximum standard deviation in experimental results of <= %0.2f;\n\t- and a resolution of <= %0.2f Angstroms.\nis %d.\n" % (MAX_NUMRES_PROTHERM, MAX_STANDARD_DEVIATION, MAX_RESOLUTION, len(experimentIDs)))
	colortext.message("The number of small-to-large mutations in the database is %d, described in %d experiments." % (len(ddGdb.execute(SQLQueries["SmallToLarge"])), len(ddGdb.execute(SQLQueries["SmallToLargeD"]))))
	colortext.message("The number of large-to-small mutations in the database is %d, described in %d experiments." % (len(ddGdb.execute(SQLQueries["LargeToSmall"])), len(ddGdb.execute(SQLQueries["LargeToSmallD"]))))
	colortext.message("The number of experiments discounted as they involved more than one chain is %d.\n" % moreThanOneChain)
	return
	raise colortext.Exception("Skip this") #todo
	for experimentID in experimentIDs:
		inserter.add(experimentID, PredictionSet, CommandName, KeepHETATMLines)


#runner = JobTestRunner()
#runner.runPrediction({dbfields.ExperimentID : 8468})
		
def main():
	#All results were produced with revision 32231 of rosetta, and revision 32257 of the rosetta database.
	jobID = 69064

	if True:
		inserter = JobInserter()
		#inserter.addAllEligibleProTherm("testrun", None, None, None)
		inserter.add(jobID, "testrun", "Kellogg:10.1002/prot.22921:protocol16:32231", False)
		sys.exit(0)
	
	#params = {
	#	"PDB_ID" : 
	#}
	
	pdbID = params["pdbID"]
	pdbfile = "%s.pdb" % pdbID
	lstfile = "%s.lst" % pdbID
	
	# Create pdb file
	F = open(os.path.join("test", pdbfile), "w")
	F.write(params[dbfields.StrippedPDB])
	F.close()

	# Create lst file
	F = open(os.path.join("test", lstfile), "w")
	F.write(pdbfile)
	F.close()
			
	mutations = params["remappedmutations"]
	pdb = PDB(params[dbfields.StrippedPDB].split("\n"))
	
	# Create resfile
	resfile = os.path.join("test", "%s.resfile" % pdbID)
	F = open(resfile, "w")
	F.write(inserter.createResfile(pdb, mutations))
	F.close()
	
	preminParams = {
		"EXECUTABLE" : "/home/oconchus/ddg/r3.3/minimize_with_cst.static.linuxgccrelease",
		"INPUT_PDB" : pdbfile,
		"INPUT_PDB_LIST" : lstfile,
		"RESFILE" : resfile,
		"DATABASE" : "/home/oconchus/ddgsrc/r3.3/rosetta_database/",
		}
	runPreminimization(preminParams, params, pdbID)
	
	ddGdb = ddgproject.ddGDatabase()
	scores = [result for result in ddGdb.callproc("GetScores", parameters = (jobID, ), cursorClass = ddgproject.StdCursor)]
	print(scores)


def createConstraintsFile(preminimizationLog):
	'''This does the work of convert_to_cst_file.sh'''
	print("Creating constraints file")
	preminimizationOutput = open(preminimizationLog, "r")
	constaints = []
	while True:
		line = preminimizationOutput.readline()
		if line == "":
			break
		elif line.startswith("c-alpha"):
			line = line.split()
			constaints.append("AtomPair CA %s CA %s HARMONIC %s %s" % (line[5], line[7], line[9], line[12]))
	preminimizationOutput.close()
	return join(constaints, "\n")

def runPreminimization(preminParams, params, pdbID):

	preminimizationLog = "test/preminimization.log"
	cstfile = "test/constraints.cst"
	
	preminimization = params[dbfields.Command]['Preminimization'] 
	args = [p % preminParams for p in preminimization]
	print("Running preminimization")
	print(join(args, " "))
	
	inputPDB = os.path.join("test", "min_cst_0.5.%s.pdb" % pdbID)
	if not(os.path.exists(inputPDB)):
		preminimizationOutput = open(preminimizationLog, "w")
		p = subprocess.Popen(args, cwd = os.path.join(os.getcwd(), "test"), stdout = preminimizationOutput, close_fds = True)
		returncode = None
		while returncode == None:
			returncode = p.poll()
			time.sleep(5)
		if returncode != 0:
			raise Exception("Preminimization failed.")
	
		os.rename(os.path.join("test", "min_cst_0.5.%s_0001.pdb" % pdbID), inputPDB)
	
		# Create the constraints file
		constraintsFile = open(cstfile, "w")
		constraintsFile.write(createConstraintsFile(preminimizationLog))
		constraintsFile.close()
	
	#"DATABASE" : "/home/oconchus/ddgsrc/r32257/rosetta_database/",

	ddGParams = {
		"EXECUTABLE" : "/home/oconchus/ddg/r32231/fix_bb_monomer_ddg.linuxgccrelease",
		"INPUT_PDB" : os.path.join(os.getcwd(), inputPDB),
		"RESFILE" : os.path.join(os.getcwd(), preminParams["RESFILE"]),
		"DATABASE" : "/home/oconchus/ddgsrc/r32257/rosetta_database/",
		"CONSTRAINTS_FILE" : os.path.join(os.getcwd(), cstfile),
		}
	
	protocols = range(0,21)
	ddgexecutable = ""
	
	commonstr = [
		'%(EXECUTABLE)s',
		'-in:file:s', '%(INPUT_PDB)s',
		'-resfile', '%(RESFILE)s',
		'-database', '%(DATABASE)s',
		'-ignore_unrecognized_res',
		'-in:file:fullatom',
		'-constraints::cst_file', '%(CONSTRAINTS_FILE)s']
	
	softrep = ['-score:weights', 'soft_rep_design']
	hardrep = "-score:weights standard –score:patch score12"
	minnohardrep = " -ddg::minimization_scorefunction standard –ddg::minimization_patch score12"
		
	protocols1617 = [
		'-ddg::weight_file', 'soft_rep_design',
		'-ddg::iterations', '1', # todo: 50
		'-ddg::local_opt_only', 'false',
		'-ddg::min_cst', 'true',
		'-ddg::mean', 'false',
		'-ddg::min', 'true',
		'-ddg::sc_min_only', 'false', # Backbone and sidechain minimization
		'-ddg::ramp_repulsive', 'true', 
		'-ddg::minimization_scorefunction', 'standard',
		'-ddg::minimization_patch', 'score12']
	
	protocols[16] = commonstr + softrep +  protocols1617
	
	args = [p % ddGParams for p in protocols[16]]
	#for i in protocols[16]:
	#	print(i % ddGParams)
	
	ddGLog = "test/ddG.log"
	ddGOutput = open(ddGLog, "w") 
	print("Running ddG")
	print(join(args, " "))
	p = subprocess.Popen(args, cwd = os.path.join(os.getcwd(), "test"), stdout = ddGOutput , close_fds = True)
	returncode = None
	while returncode == None:
		returncode = p.poll()
		time.sleep(5)
	if returncode != 0:
		raise Exception("ddG failed with return code %d" % returncode)

	print(parseResults("test/ddG.log", "test/ddg_predictions.out"))
	return

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

# Misc functions
def checkForInsertionCodeEntries():
	'''This function identifies mutations which were not included in the initial data import
	IDs are (5438, 5439, 5440, 5441, 8060, 13083, 13084, 14584, 15616, 15617, 15626, 15627, 15628, 15629, 15638, 15639, 15640, 15641)
	Nnoe of these have ddG values but we may need to pay attention to these if the ProTherm database is updated.'''
	
	idmatcher = re.compile("NO.\s+(\d+)")
	mutationmatcher		= re.compile("(\w\s+\d+\s+\w[,]*)+")
	newmutationmatcher	= re.compile("(\w\s+\d+\w?\s+\w[,]*)+")
	F = open("rawdata/ProTherm.dat")
	id = None
	ids = []
	print_ddG= False
	for line in F:
		if line.startswith("NO."):
			print_ddG= False
			id = idmatcher.match(line).groups(1)
		if line.startswith("MUTATION"):
			data = line[8:].strip()
			if data and data != "wild**" and data != "wild*" and data != "wild":
				matches = mutationmatcher.match(data)
				if not(matches):
					print_ddG = True
					assert(len(id) == 1)
					ids.append(int(id[0]))
					matches = newmutationmatcher.match(data)
					if not(matches):
						colortext.printf("%s: %s" % (id, line.strip()), "orange")
					else:
						colortext.printf("%s: %s" % (id, line.strip()), "yellow")
		if print_ddG and line.startswith("ddG") and not line.startswith("ddG_"):
			ddG = line[3:].strip()
			if ddG:
				colortext.printf("\tddG: '%s'" % ddG, "green")
			else:
				colortext.printf("\tNo ddG:", "red")
		
	F.close()
	print(tuple(ids))
	
main()
	