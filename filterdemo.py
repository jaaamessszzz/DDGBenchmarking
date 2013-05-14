import sys
import time

sys.path.insert(0, "..")
#sys.path.insert(0, "ddglib")
from tools import colortext
from tools.deprecated import rosettadb
from tools.debug.profile import ProfileTimer
from ddglib import dbapi, ddgdbapi, help
from ddglib.ddgfilters import *

#help.help()

def simpleRunExample(self):
	# Step 1: Open a database connection
	ddGdb = ddgdbapi.ddGDatabase()

	# Step 2: Select database records
	sr = StructureResultSet(ddGdb, AdditionalIDs = ['2BQC', '1LAW', '1LHH', '1LHI'])
	
	# Step 3: Add filters
	sr.addFilter(StructureFilter.TotalBFactors(0,16) | StructureFilter.WithNullResolution(True))
	
	# Step 4: Retrieve full database records. 
	# results will be a list each of whose elements is a dict representing a database record.
	results = sr.getFilteredResults()
	
	# Step 5: Optionally print summary
	print("\nSummary: %s\n" % sr)


# Create a result set
class Examples:
	
	@staticmethod
	def printOutput(resultset):
		print("Applying filters")
		results = resultset.getFilteredResults()
		print("After application")
		print("\nSummary: %s\n" % resultset)

	@staticmethod
	def openDB():
		if not globals().get("ddGdb"):
			globals()["ddGdb"] = ddgdbapi.ddGDatabase()

	@staticmethod
	def help():
		help.ShowDatabaseStructure()
		help.ShowResultSet()
		help.ShowFilter()

	# UnionFilter examples

	@staticmethod
	def unionFilterExample1():
		print("** All structures with null OR non-null resolution **")
		pt = ProfileTimer()
		pt.start("Open DB")
		Examples.openDB()
		pt.start("Create StructureResultSet")
		sr = StructureResultSet(ddGdb)
		pt.start("addFilter")
		sr.addFilter(StructureFilter.WithNullResolution(False) | StructureFilter.WithNullResolution(True))
		pt.start("printOutput")
		Examples.printOutput(sr)
		pt.stop()
		print("\nProfile:\n\n%s\n" % pt.getProfile(html = False))

	@staticmethod
	def unionFilterExample2():
		print("** All structures with null AND non-null resolution**") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.WithNullResolution(False))
		sr.addFilter(StructureFilter.WithNullResolution(True))
		Examples.printOutput(sr)


	# StructureResultSet examples

	@staticmethod
	def allStructures():
		'''Select all Structure records.'''
		print("** All strucures **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb)
		Examples.printOutput(sr)

	@staticmethod
	def getStructuresWithNullResolutionSQL():
		print("** All structures with null resolution **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb, SQL = "WHERE Resolution IS NULL")
		Examples.printOutput(sr)

	@staticmethod
	def getStructuresWithNullResolutionFilter():
		print("** All structures with null resolution **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.WithNullResolution(True))
		Examples.printOutput(sr)

	@staticmethod
	def pickSpecific():
		'''Select four specific Structure records and apply a filter.''' 
		print("** 4 specific structures **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb, AdditionalIDs = ['2BQC', '1LAW', '1LHH', '1LHI'])
		sr.addFilter(StructureFilter.TotalBFactors(0,16) | StructureFilter.WithNullResolution(True))
		Examples.printOutput(sr)

	@staticmethod
	def getStructuresInResolutionRange():
		print("** All structures with null resolution **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.Resolution(1, 2))
		Examples.printOutput(sr)
	
	@staticmethod
	def getStructuresWithUniProtIDs():
		print("** All structures with null resolution **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.WithUniProtIDs(["P0A7Y4"], ["RNH_ECOLI", "RNP30_RANPI"]))
		Examples.printOutput(sr)

	@staticmethod
	def getStructuresFilteredByStructures():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **") 
		Examples.openDB()
		
		sr1 = StructureResultSet(ddGdb, SQL = "WHERE PDB_ID LIKE %s", parameters = "1A%")
		Examples.printOutput(sr1)
		
		sr2 = StructureResultSet(ddGdb, SQL = "WHERE PDB_ID LIKE %s", parameters = "1AY%")
		Examples.printOutput(sr2)
		
		sr = sr1.filterBySet(sr2)
		Examples.printOutput(sr)


	# ExperimentResultSet examples

	@staticmethod
	def getExperimentsWithSQL():
		'''Select all Structure records.'''
		print("** All structures **") 
		Examples.openDB()
		er = ExperimentResultSet(ddGdb, SQL = "WHERE Structure LIKE %s", parameters = "1A%")
		Examples.printOutput(er)
		print(er.structure_map.keys())
		
		er.addFilter(StructureFilter.Resolution(1, 1.7))
		
		Examples.printOutput(er)
		
	@staticmethod
	def getExperimentsFilteredByStructures():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **") 
		Examples.openDB()
		
		sr = StructureResultSet(ddGdb, SQL = "WHERE PDB_ID LIKE %s", parameters = "1AY%")
		Examples.printOutput(sr)
		
		er = ExperimentResultSet(ddGdb, SQL = "WHERE Structure LIKE %s", parameters = "1A%")
		Examples.printOutput(er)
		
		er = er.filterBySet(sr)
		Examples.printOutput(er)
		
		er = ExperimentResultSet(ddGdb, SQL = "WHERE Structure LIKE %s", parameters = "1AY%")
		Examples.printOutput(er)
		
		#print(er.structure_map.keys())
		
		er.addFilter(StructureFilter.Resolution(1, 1.7))
		
		Examples.printOutput(er)
	
	@staticmethod
	def getExperimentsFilteredBySource():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **") 
		Examples.openDB()
		
		er = ExperimentResultSet(ddGdb)
		Examples.printOutput(er)
		
		er.addFilter(ExperimentFilter.OnSource(ExperimentFilter.ProTherm))
		
		Examples.printOutput(er)

	
	@staticmethod
	def getExperimentsFilteredByMutationSize():
		'''Select all Structure records.'''
		print("** Experiments filtered by mutation size **") 
		Examples.openDB()
		
		er = ExperimentResultSet(ddGdb)
		Examples.printOutput(er)
		
		er.addFilter(ExperimentFilter.MutationsBetweenAminoAcidSizes(ExperimentFilter.large, ExperimentFilter.large))
		
		Examples.printOutput(er)

	@staticmethod
	def getExperimentsFilteredByAminoAcids1():
		'''Select all Structure records.'''
		print("** Experiments filtered by residue (from ALA) **") 
		Examples.openDB()
		
		er = ExperimentResultSet(ddGdb)
		Examples.printOutput(er)
		
		er.addFilter(ExperimentFilter.MutationsBetweenAminoAcids('ALA', 'G'))
		
		Examples.printOutput(er)

	@staticmethod
	def getExperimentsFilteredByAminoAcids2():
		'''Select all Structure records.'''
		print("** Experiments filtered by residue (from ALA) **") 
		Examples.openDB()
		
		er = ExperimentResultSet(ddGdb)
		Examples.printOutput(er)
		
		er.addFilter(ExperimentFilter.MutationsBetweenAminoAcids('A', 'GLY'))
		
		Examples.printOutput(er)

	@staticmethod
	def getExperimentsFilteredBySourceAndResolution():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **") 
		Examples.openDB()
		
		er = ExperimentResultSet(ddGdb)
		Examples.printOutput(er)
		
		er.addFilter(ExperimentFilter.OnSource(ExperimentFilter.ProTherm))
		
		Examples.printOutput(er)
		
		er.addFilter(StructureFilter.Resolution(1, 2))
		Examples.printOutput(er)
		
		
	# PredictionResultSet examples

	@staticmethod
	def getPredictionsWithSQL():
		'''Select all Structure records.'''
		print("** All structures **") 
		Examples.openDB()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s AND ID=12595", parameters = "testrun")
		Examples.printOutput(pr)

	@staticmethod
	def getPredictionsUsingMultipleFilters():
		'''This demonstrates the use of multiple filters.'''
		print("** Multiple filter example **") 
		Examples.openDB()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s", parameters = "testrun")
		pr.addFilter(StructureFilter.Techniques(StructureFilter.XRay))
		pr.addFilter(StructureFilter.Resolution(1, 1.5) | StructureFilter.Resolution(3.9, 4))
		pr.addFilter(StructureFilter.TotalBFactors(0, 10))
		Examples.printOutput(pr)

	@staticmethod
	def getPredictionsUsingMultipleFilters_Speed():
		'''This demonstrates how slow separate filters are.'''
		print("** Multiple filter example **") 
		Examples.openDB()

		import time
		
		t1 = time.time()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s", parameters = "testrun")
		pr.addFilter(StructureFilter.Techniques(StructureFilter.XRay))
		pr.addFilter(StructureFilter.Resolution(1, 1.5) | StructureFilter.Resolution(3.9, 4))
		pr.addFilter(StructureFilter.TotalBFactors(0, 10))
		Examples.printOutput(pr)
		t2 = time.time()

		print("Time taken: %0.2fs" % (t2 - t1))

		t1 = time.time()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s", parameters = "testrun")
		sf = StructureFilter()
		sf.setTechniques(StructureFilter.XRay)
		sf.setResolution(1, 1.5)
		sf.setTotalBFactors(0, 10)
		pr.addFilter(sf | StructureFilter.Resolution(3.9, 4))
		Examples.printOutput(pr)
		t2 = time.time()
		print("Time taken: %0.2fs" % (t2 - t1))

	@staticmethod
	def getPredictionsUsingMultipleFilters2():
		print("** Multiple filter example **") 
		Examples.openDB()
		
	@staticmethod
	def showResultSetOperations():
		'''Demonstrates how to union, intersect, subtract, and XOR ResultSets.'''
		print("\n** ResultSet SR1 **\n")
		Examples.openDB()
		sr1 = StructureResultSet(ddGdb)
		sr1.addFilter(StructureFilter.Resolution(1, 1.3))
		Examples.printOutput(sr1)
		
		print("\n** ResultSet SR2 **\n")
		sr2 = StructureResultSet(ddGdb)
		sr2.addFilter(StructureFilter.Resolution(2, 2.3))
		Examples.printOutput(sr2)
		
		print("\n** ResultSet SR3 **\n")
		sr3 = StructureResultSet(ddGdb)
		sr3.addFilter(StructureFilter.Resolution(1.2, 2))
		Examples.printOutput(sr3)
		
		print("\n** ResultSet union - SR1 | SR2 **\n")
		srUnion = sr1 | sr2
		print(join(srUnion._log, "\n"))
		
		print("\n** ResultSet union - SR1 - SR3 **\n")
		
		srUnion = sr1 | sr3
		print(join(srUnion._log, "\n"))

		print("\n** ResultSet intersection - SR1 & SR3 **\n")
		
		srIntersection = sr1 & sr3
		print(join(srIntersection._log, "\n"))
		
		print("\n** ResultSet intersection sanity check **\n")
		
		sr4 = StructureResultSet(ddGdb)
		sr4.addFilter(StructureFilter.Resolution(1.2, 1.3))
		Examples.printOutput(sr4)
		
		print("\n** ResultSet difference - SR1 - SR3 **\n")
		
		#srDifference = sr1 - sr3
		srDifference = sr1 / sr3
		print(join(srDifference._log, "\n"))

		print("\n** ResultSet exclusive or - SR1 ^ SR3 **\n")
		
		srXOR = sr1 ^ sr3
		print(join(srXOR._log, "\n"))
	
	@staticmethod
	def showAllEligibleProTherm(PredictionSet, ProtocolID, KeepHETATMLines):
		#inserter = JobInserter()
		colortext.printf("\nAdding ProTherm mutations to %s prediction set." % PredictionSet, "lightgreen")
		#ddGdb = ddgdbapi.ddGDatabase()
		
		MAX_RESOLUTION = 2.1
		MAX_NUMRES_PROTHERM = 350
		MAX_STANDARD_DEVIATION = 1.0

		Examples.openDB()
		
		import time
		if False:
			t1 = time.time()
			er1 = ExperimentResultSet(ddGdb)
			er1.addFilter(ExperimentFilter.OnSource(ExperimentFilter.ProTherm))
			er1.addFilter(ExperimentFilter.NumberOfMutations(1, 1))
			er1.addFilter(ExperimentFilter.NumberOfChains(1, 1))
			er1.addFilter(ExperimentFilter.StandardDeviation(None, MAX_STANDARD_DEVIATION))
			er1.addFilter(StructureFilter.Resolution(None, MAX_RESOLUTION))
			er1.addFilter(StructureFilter.Techniques(StructureFilter.XRay))
			Examples.printOutput(er1)
			t2 = time.time()
			print(t2 - t1)
		
		# This method usually takes around 65% of the time as the method above 
		t1 = time.time()
		ef1 = ExperimentFilter()
		ef1.setSource(ExperimentFilter.ProTherm)
		er1 = ExperimentResultSet(ddGdb)
		er1.addFilter(ExperimentFilter.OnSource(ExperimentFilter.ProTherm))
		Examples.printOutput(er1)
		ef1.setNumberOfMutations(1, 1)
		ef1.setNumberOfChains(1, 1)
		ef1.setStandardDeviation(None, MAX_STANDARD_DEVIATION)
		sf1 = StructureFilter()
		sf1.setResolution(None, MAX_RESOLUTION)
		sf1.setTechniques(StructureFilter.XRay)
		er1 = ExperimentResultSet(ddGdb)
		er1.addFilter(ef1)
		er1.addFilter(sf1)
		Examples.printOutput(er1)
		t2 = time.time()
		print(t2 - t1)
		
		experimentIDs = sorted(list(er1.getFilteredIDs()))
		colortext.message("\nThe number of unique ProTherm experiments with:\n\t- one mutation;\n\t- structures solved by X-ray diffraction and with <= %d residues;\n\t- a maximum standard deviation in experimental results of <= %0.2f;\n\t- and a resolution of <= %0.2f Angstroms.\nis %d.\n" % (MAX_NUMRES_PROTHERM, MAX_STANDARD_DEVIATION, MAX_RESOLUTION, len(experimentIDs)))
		ddG_connection = dbapi.ddG()
		count = 0
		sys.exit(0)
		print("")
		for experimentID in experimentIDs:
			ddG_connection.addPrediction(experimentID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True)
			count += 1
			if count >= 10:
				colortext.write(".")
				colortext.flush()
				count = 0
		print("")
		
	@staticmethod
	def testAnalysis():
		ddG_connection = dbapi.ddG()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet='kellogg16-A' AND Status='done' LIMIT 2000")
		ddG_connection.analyze(pr)
	
	@staticmethod
	def testAnalysis2():
		ddG_connection = dbapi.ddG()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet='lizsettest1' AND Status='done' LIMIT 2000")
		ddG_connection.analyze(pr)

		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet='lizsettest1' AND Status='done' LIMIT 2000")
		pr.addFilter(ExperimentFilter.MutationsBetweenAminoAcidSizes(ExperimentFilter.large, ExperimentFilter.small))
		Examples.printOutput(pr)
		ddG_connection.analyze(pr)

		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet='lizsettest1' AND Status='done' LIMIT 2000")
		pr.addFilter(ExperimentFilter.MutationsBetweenAminoAcidSizes(ExperimentFilter.small, ExperimentFilter.large))
		Examples.printOutput(pr)
		ddG_connection.analyze(pr)

	@staticmethod
	def testPublications():
		ddG_connection = dbapi.ddG()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE ID >= 12804 and ID <= 12903")
		er = ExperimentResultSet(ddGdb, SQL = "WHERE ID >= 73534 and ID <= 73561")
		ddG_connection.getPublications(pr)
		ddG_connection.getPublications(er)
	
	@staticmethod
	def addAllMutationsForAGivenPDB1():
		import common.pdb
		import common.rosettahelper
		
		ddG_connection = dbapi.ddG()
		#ddG_connection.addPDBtoDatabase("3K0NA_lin.pdb")
		opdb = common.pdb.PDB("3K0NA_lin.pdb")
		common.rosettahelper.ROSETTAWEB_SK_AAinv
		count = 1
		for chainresidueid, wt in sorted(opdb.ProperResidueIDToAAMap().iteritems()):
			chain = chainresidueid[0]
			residueid = chainresidueid[1:].strip()
			allotherAAs = sorted([aa for aa in common.rosettahelper.ROSETTAWEB_SK_AAinv.keys() if aa != wt])
			for otherAA in allotherAAs: 
				ms = dbapi.MutationSet()
				ms.addMutation(chain, residueid, wt, otherAA)
				print("3K0NA_lin", ms, ms.getChains(), count, 0)
				ddG_connection.createDummyExperiment("3K0NA_lin", ms, ms.getChains(), count, 0, ExperimentSetName = "DummySource")
				count += 1

	@staticmethod
	def addAllMutationsForAGivenPDB2():
		import common.pdb
		import common.rosettahelper
		
		ddG_connection = dbapi.ddG()
		#ddG_connection.addPDBtoDatabase("3K0On_lin.pdb")
		opdb = common.pdb.PDB("3K0On_lin.pdb")
		common.rosettahelper.ROSETTAWEB_SK_AAinv
		count = 3098
		for chainresidueid, wt in sorted(opdb.ProperResidueIDToAAMap().iteritems()):
			chain = chainresidueid[0]
			residueid = chainresidueid[1:].strip()
			allotherAAs = sorted([aa for aa in common.rosettahelper.ROSETTAWEB_SK_AAinv.keys() if aa != wt])
			for otherAA in allotherAAs: 
				ms = dbapi.MutationSet()
				ms.addMutation(chain, residueid, wt, otherAA)
				print("3K0On_lin", ms, ms.getChains(), count, 0)
				ddG_connection.createDummyExperiment("3K0On_lin", ms, ms.getChains(), count, 0, ExperimentSetName = "DummySource")
				count += 1
				
	@staticmethod
	def addAllMutationsForAGivenPDB3():
		import common.pdb
		import common.rosettahelper
		
		Examples.openDB()
		ddG_connection = dbapi.ddG()
		#ddG_connection.addPDBtoDatabase("pdbs/3K0NB_lin.pdb", UniProtAC = "P62937", UniProtID = "PPIA_HUMAN")
		opdb = common.pdb.PDB("pdbs/3K0NB_lin.pdb")
		common.rosettahelper.ROSETTAWEB_SK_AAinv
		
		results = ddGdb.execute('''SELECT SourceID FROM ExperimentScore INNER JOIN Experiment ON ExperimentScore.ExperimentID = Experiment.ID WHERE Source="DummySource"''', cursorClass=ddgdbapi.StdCursor) 
		assert(results)
		highestID = max([int(r[0]) for r in results])
		count = highestID + 1

		for chainresidueid, wt in sorted(opdb.ProperResidueIDToAAMap().iteritems()):
			chain = chainresidueid[0]
			residueid = chainresidueid[1:].strip()
			allotherAAs = sorted([aa for aa in common.rosettahelper.ROSETTAWEB_SK_AAinv.keys() if aa != wt])
			for otherAA in allotherAAs: 
				ms = dbapi.MutationSet()
				ms.addMutation(chain, residueid, wt, otherAA)
				print("3K0NB_lin", ms, ms.getChains(), count, 0, chain, wt, residueid, otherAA)
				ddG_connection.createDummyExperiment("3K0NB_lin", ms, ms.getChains(), count, 0, ExperimentSetName = "DummySource")
				count += 1
				
	@staticmethod
	def addLinsJobs(PredictionSet, ProtocolID):
		colortext.printf("\nAdding Lin's mutations to %s prediction set." % PredictionSet, "lightgreen")
		KeepHETATMLines = False
		Examples.openDB()
		
		# Filter by the DummySource set of experiments
		er1 = ExperimentResultSet(ddGdb)
		ef1 = ExperimentFilter()
		ef1.setSource(ExperimentFilter.DummySource)
		er1.addFilter(ef1)
		
		# Filter by the particular PDB
		sr = StructureResultSet(ddGdb, 'WHERE PDB_ID="3K0NB_lin"')
		er1 = ExperimentResultSet.fromIDs(ddGdb, er1.getFilteredIDs()).filterBySet(sr)
		Examples.printOutput(er1)
		
		experimentIDs = sorted(list(er1.getFilteredIDs()))
		colortext.message("\nThe number of unique experiments is %d.\n" % len(experimentIDs))
		ddG_connection = dbapi.ddG()
		count = 0
		for experimentID in experimentIDs:
			ddG_connection.addPrediction(experimentID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True)
			count += 1
			if count >= 10:
				colortext.write(".")
				colortext.flush()
				count = 0
		print("")

	@staticmethod
	def runLizsSet(PredictionSet, ProtocolID):
		colortext.printf("\nAdding Liz's data set to %s prediction set." % PredictionSet, "lightgreen")
		KeepHETATMLines = False
		Examples.openDB()
		
		# Filter by the DummySource set of experiments
		er1 = ExperimentResultSet(ddGdb)
		ef1 = ExperimentFilter()
		ef1.setSource(ExperimentFilter.LizKellogg)
		er1.addFilter(ef1)
		Examples.printOutput(er1)
		
		experimentIDs = sorted(list(er1.getFilteredIDs()))
		colortext.message("\nThe number of unique experiments is %d.\n" % len(experimentIDs))
		ddG_connection = dbapi.ddG()
		count = 0
		for experimentID in experimentIDs:
			ddG_connection.addPrediction(experimentID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True)
			count += 1
			if count >= 10:
				colortext.write(".")
				colortext.flush()
				count = 0
		print("")

if False:
	import analysis
	analyzer = analysis.Analyzer("AllExperimentsProtocol16")
	analyzer.AddPublishedDDGsToAnalysisTables()
	analyzer.plot(analysis.Analyzer.correlation_coefficient, "Kellogg.rr", table_names = ["Kellogg"])
	#			"kellogg.txt", "Kellogg")
	for table_name, a_table in sorted(analyzer.analysis_tables.iteritems()):
		print(a_table)
		print(table_name)
	
#print(analysis.AnalysisPoint.headers)
#print(analysis_tables)
#print(analysis.AnalysisPoint.headers)
#print(analysis_tables["Kellogg"])

#ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "AllExperimentsProtocol16", "Kellogg:10.1002/prot.22921:protocol16:32231", False, StoreOutput = False, Description = {}, InputFiles = {}, testonly = False)
	
#ddG_connection = dbapi.ddG()
#ddG_connection.addPDBtoDatabase(pdbID = "1FKJ")

ddG_connection = dbapi.ddG()

Examples.help()
#ddG_connection = dbapi.ddG()
#ddG_connection.dumpData("testzip-13103.zip", 13103)
#Examples.addAllMutationsForAGivenPDB3()
Examples.unionFilterExample1()
#Examples.showAllEligibleProTherm("test", "test", False)
#Examples.addLinsJobs("lin-3K0NB", "Kellogg:10.1002/prot.22921:protocol16:32231")
#Examples.testAnalysis2()
#Examples.runLizsSet("lizsettest1", "Kellogg:10.1002/prot.22921:protocol16:32231")
#Examples.testAnalysis()
#Examples.testPublications()

#Examples.getExperimentsFilteredByStructures()
#Examples.getStructuresFilteredByStructures()
#Examples.showResultSetOperations()
#Examples.getStructuresWithUniProtIDs()
#Examples.getExperimentsFilteredByStructures()
#Examples.getExperimentsFilteredBySource()
#Examples.getExperimentsFilteredBySourceAndResolution()
#Examples.getExperimentsFilteredByMutationSize()
#Examples.getExperimentsFilteredByAminoAcids1()
#Examples.getExperimentsFilteredByAminoAcids2()
#Examples.help()
#Examples.getPredictionsUsingMultipleFilters_Speed()
#help.ShowFilter()
#Examples.allStructures()
#Examples.pickSpecific()
