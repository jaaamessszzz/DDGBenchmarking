import common.ddgproject
from ddglib.help import show_tables
from ddglib.ddgfilters import *


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
	def allStructures():
		'''Select all Structure records.'''
		print("** All strucures **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb)
		Examples.printOutput(sr)

	@staticmethod
	def pickSpecific():
		'''Select four specific Structure records and apply a filter.''' 
		print("** 4 specific structures **") 
		Examples.openDB()
		sr = StructureResultSet(ddGdb, AdditionalIDs = ['2BQC', '1LAW', '1LHH', '1LHI'])
		sr.addFilter(StructureFilter.TotalBFactors(0,16) | StructureFilter.WithNullResolution(True))
		Examples.printOutput(sr)

show_tables()
if False:
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
import sys
sys.exit(0)

#ddGdb = common.ddgproject.ddGDatabase()
Examples.allStructures()
Examples.pickSpecific()
