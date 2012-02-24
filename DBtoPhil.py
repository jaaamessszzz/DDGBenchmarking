import sys
import os
from string import join
import math
import pickle
import re
sys.path.insert(0, "common")
import ddgproject
import colortext

stdCursor = ddgproject.StdCursor

db = ddgproject.ddGDatabase()
experiments = db.execute("SELECT * FROM Experiment")

print(40 * "*")
step = len(experiments) / 40
count = 0

headers = [
	"KortemmeLabID", "Source", "PDB ID", "UniProt AC/IDs", "Resolution (A)",
	"Experimental Techniques", "BFactor average", "BFactor standard deviation",
	"Mutant structures", "Interfaces", "Chains", "Mutations", "Scores",
]

F = open("KortemmeLabInterfaces-2012-01-13.txt", "w")
F.write("%s\n" % join(headers, "\t"))
for experiment in experiments:
	ID = experiment["ID"] 
	SOURCE = experiment["Source"] 
	PDB_ID = experiment["Structure"] 
	
	uniprot = sorted(db.execute("SELECT UniProtKBMapping.UniProtKB_AC AS AC, UniProtKB.UniProtKB_ID AS ID FROM UniProtKBMapping INNER JOIN UniProtKB WHERE UniProtKBMapping.UniProtKB_AC=UniProtKB.UniProtKB_AC AND PDB_ID=%s", parameters = (PDB_ID,), cursorClass = stdCursor))
	UNIPROT = ["(%s, %s)" % mapping for mapping in uniprot]
	UNIPROT = join(UNIPROT, "")
	
	pdbdetails = db.execute("SELECT Resolution, Techniques, BFactors FROM Structure WHERE PDB_ID=%s", parameters = (PDB_ID,))[0]
	RESOLUTION = pdbdetails["Resolution"] or "N/A"
	TECHNIQUES = (pdbdetails["Techniques"]).replace(";", ",")
	BFactors = pickle.loads(pdbdetails["BFactors"])["Total"]
	BFACTOR_AVG = BFactors[0] 
	BFACTOR_STDDEV = BFactors[1] 
	
	mutants = sorted(db.callproc("GetMutants", ID, cursorClass = stdCursor))
	for m in mutants:
		assert(len(m) == 1)
	MUTANTS = [m[0] for m in mutants]
	MUTANTS = join(MUTANTS, ",")
	
	interfaces = db.callproc("GetInterfaces", ID, cursorClass = stdCursor)
	if len(interfaces) > 1:
		print("INTERFACES", interfaces)
		sys.exit(0)
	for i in interfaces:
		assert(len(i) == 1)
	INTERFACES = [i[0] for i in interfaces]
	INTERFACES = join(INTERFACES, ",")
	
	chains = db.callproc("GetChains", ID, cursorClass = stdCursor)
	for c in chains:
		assert(len(c) == 1)
	CHAINS = [c[0] for c in chains]
	CHAINS = join(CHAINS, ",")
	
	mutations = db.callproc("GetMutations", ID, cursorClass = stdCursor)
	MUTATIONS = []
	for m in mutations:
		MUTATIONS.append("(%s, %s, %s, %s)" % m)
	MUTATIONS = join(MUTATIONS, ",")
	
	scores = db.callproc("GetScores", ID)
	SCORES = []
	for s in scores:
		SCORES.append("(%s, %s, %s)" % (s["SourceID"], s["ddG"], s["NumberOfMeasurements"]))
	SCORES = join(SCORES, ",")
	
	record = map(str, [
		ID, SOURCE, PDB_ID, UNIPROT, RESOLUTION, 
		TECHNIQUES, BFACTOR_AVG, BFACTOR_STDDEV,
		MUTANTS, INTERFACES, CHAINS, MUTATIONS, SCORES])
	assert(len(headers) == len(record)) 
	record = join(record, "\t")
	F.write("%s\n" % record)
	
	count += 1
	if count >= step:
		sys.stdout.write(".")
		sys.stdout.flush()
		count = 0

sys.stdout.write("\n")
F.close()
	