import common.ddgproject
from ddglib.ddgfilters import *

def printOutput(resultset):
	print("Applying filters")
	results = resultset.getFilteredResults()
	print("After application")
	print("\nSummary: %s\n" % resultset)
	for r in results:
		print(r["PDB_ID"], pickle.loads(r["BFactors"])["Total"])

# Open a database connection
ddGdb = common.ddgproject.ddGDatabase()

# Create a result set
sr = StructureResultSet(ddGdb, AdditionalIDs = ['2BQC', '1LAW', '1LHH', '1LHI'])
sr = StructureResultSet(ddGdb)
sr.addFilter(StructureFilter.TotalBFactors(0,16) | StructureFilter.WithNullResolution(True))

printOutput(sr)

