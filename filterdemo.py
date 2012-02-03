import common.ddgproject
from ddglib import help
from ddglib.ddgfilters import *

help.help()

def simpleRunExample(self):
	# Step 1: Open a database connection
	ddGdb = common.ddgproject.ddGDatabase()

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
			globals()["ddGdb"] = common.ddgproject.ddGDatabase()

	@staticmethod
	def help():
		pass
		#help.ShowDatabaseStructure()
		#help.ShowResultSet()
		#help.ShowFilter()

	# UnionFilter examples

	@staticmethod
	def unionFilterExample1():
		print("** All structures with null OR non-null resolution **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.WithNullResolution(False) | StructureFilter.WithNullResolution(True))
		Examples.printOutput(sr)

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
		#Examples.printOutput(pr)
		pr.addFilter(StructureFilter.Resolution(1, 1.5) | StructureFilter.Resolution(3.9, 4))
		#Examples.printOutput(pr)
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
		
#ddGdb = common.ddgproject.ddGDatabase()
Examples.help()
Examples.getPredictionsUsingMultipleFilters_Speed()
help.ShowFilter()
#Examples.allStructures()
#Examples.pickSpecific()
