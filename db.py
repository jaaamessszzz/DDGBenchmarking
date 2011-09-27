#!/usr/bin/python2.4
# -*- coding: iso-8859-15 -*-

import sys, os
import MySQLdb
import MySQLdb.cursors
import traceback
import pickle
import time
from datetime import datetime, date
from string import join
import math
from httplib import HTTPConnection

sqrt = math.sqrt
DictCursor = MySQLdb.cursors.DictCursor

PDBToUniProt = {}
UniProtKBACToPDB = {}
uniprotmapping = os.path.join("rawdata", "uniprotmapping.csv")
UniqueIDs = {}
	
def readUniProtMap():
	if os.path.exists(uniprotmapping):
		F = open(uniprotmapping)
		contents = F.read().split("\n")
		F.close()
		for line in contents[1:]:
			if line:
				line = line.split("\t")
				pdbID = line[0]
				UPAC = line[1]
				UPID = line[2]
				if PDBToUniProt.get(pdbID):
					PDBToUniProt[pdbID].append((UPAC, UPID))
				else:
					PDBToUniProt[pdbID] = [(UPAC, UPID)]
				if UniProtKBACToPDB.get(UPAC):
					UniProtKBACToPDB[UPAC].append(pdbID)
				else:
					UniProtKBACToPDB[UPAC] = [pdbID]
		
class FieldNames(object):
	'''Define database fieldnames here so we can change them in one place if need be.'''
	
	PDB_ID = "PDB_ID"
	Content = "Content"
	Resolution = "Resolution"
	Source = "Source"
	Protein = "Protein"
	
	Structure = "Structure"
	Mutant = "Mutant"
	ScoreVariance = "ScoreVariance"
			
	ExperimentID = "ExperimentID"
	Chain = "Chain"
	ResidueID = "ResidueID"
	WildTypeAA = "WildTypeAA"
	MutantAA = "MutantAA"
	
	SourceID = "SourceID"
	ddG = "ddG"
	NumberOfMeasurements = "NumberOfMeasurements"
	
	Name = "Name"
	Version = "Version"
	SVNRevision = "SVNRevision"
	
	Type = "Type"
	Command = "Command"
	
	PredictionSet = "PredictionSet"
	ToolID = "ToolID"
	CommandID = "CommandID"
	KeptHETATMLines = "KeptHETATMLines"
	StrippedPDB = "StrippedPDB"
	InputFiles = "InputFiles"
	Description = "Description"
	Date = "Date"

class DBObject(object):
	dict = {}
	dbID = None
	
	def __init__(self, Description):
		pass
	
	def __getitem__(self, key):
		return self.dict[key]

	def getDBID(self):
		if not self.dbID:
			raise Exception("Cannot get the database ID of an uncommitted experiment.")
		else:
			return self.dbID
		
	def commit(self, db):
		raise Exception("Concrete function commit needs to be defined.")
	
	def __repr__(self):
		raise Exception("Concrete function  __repr__  needs to be defined.")
	
class PDBStructure(DBObject):
	
	def __init__(self, pdbID, content = None, protein = None, source = None):
		self.dict = {
			FieldNames.PDB_ID : pdbID,
			FieldNames.Content : content,
			FieldNames.Protein : protein,
			FieldNames.Source : source,
			FieldNames.Resolution : None,
		}
		
	
	def getPDBContents(self):
		d = self.dict
		id = d[FieldNames.PDB_ID]
		if len(id) != 4:
			print(id)
		assert(len(id) == 4)
		filename = os.path.join("pdbs", id + ".pdb")
		contents = None
		
		if not os.path.exists(filename):
			print("Missing file for %s %s" % (id, filename))
			c = HTTPConnection("www.rcsb.org")
			c.request("GET", "/pdb/files/%s.pdb" % id)
			response = c.getresponse()
			contents = response.read()
			c.close()
			if contents[0:6] == "<html>":
				print(contents)
				raise Exception("Error retrieving %s." % filename)				
			F = open(filename, "w")
			F.write(contents)
			F.close()		
			if not os.path.exists(filename):
				raise Exception("Error retrieving %s." % filename)
		else:
			F = open(filename, "r")
			contents = F.read()
			F.close()		
		
		resolution = None
		for line in contents.split("\n"):
			if line.startswith("REMARK   2 RESOLUTION."):
				if line[23:38] == "NOT APPLICABLE.":
					resolution = "N/A"
				elif line[31:].startswith("ANGSTROMS."):
					resolution = float(line[22:30])
				else:
					print(line)
				break

		if not resolution:
			raise Exception("Could not determine resolution for %s." % filename)
		if resolution == "N/A":
			resolution = None
		
		UniqueIDs[id] = True
			
		if not PDBToUniProt:
			readUniProtMap()
		if not PDBToUniProt.get(id):
			raise Exception("Could not find a UniProt mapping for %s in %s." % (id, uniprotmapping))
		d[FieldNames.Content] = contents
		d[FieldNames.Resolution] = resolution
			
	def commit(self, db):
		d = self.dict
		
		self.getPDBContents()
		
		results = db.execQuery("SELECT * FROM Structure WHERE PDB_ID=%s", parameters = (d[FieldNames.PDB_ID]), cursorClass = DictCursor)
		
		if results:
			assert(len(results) == 1)
			result = results[0]
			pdbID = results[0][FieldNames.PDB_ID]
			for k, v in d.iteritems():
				if k != FieldNames.PDB_ID:
					if d[k] and not(result[k]):
						SQL = "UPDATE Structure SET %s" % k
						SQL += "=%s WHERE PDB_ID=%s" 
						results = db.execQuery(SQL, parameters = (v, pdbID))
		else:
			SQL = 'INSERT INTO Structure (PDB_ID, Content, Resolution, Protein, Source) VALUES (%s, %s, %s, %s, %s);'
			vals = (d[FieldNames.PDB_ID], d[FieldNames.Content], d[FieldNames.Resolution], d[FieldNames.Protein], d[FieldNames.Source]) 
			db.execQuery(SQL, parameters = vals)			
		
	def __repr__(self):
		d = self.dict
		str = []
		str.append("%s: %s" % (FieldNames.PDB_ID, d[FieldNames.PDB_ID]))
		str.append("%s: %s" % (FieldNames.Description, d[FieldNames.Description]))
		return join(str, "\n")
	
class ExperimentSet(DBObject):
	
	def __init__(self, pdbid, source, interface = None):
		self.dict = {
			FieldNames.Structure	: pdbid,
			FieldNames.Source		: source,
			FieldNames.ScoreVariance: None,
			"Interface"				: interface,
			"Mutants"				: {},
			"Mutations"				: [],
			"ExperimentChains"		: [],
			"ExperimentalScores"	: [],
			"StdDeviation"			: None,
			"WithinStdDeviation"	: None
		}
		self.dbID = None
	
	def addMutant(self, mutant):
		self.dict["Mutants"][mutant] = True
		
	def addMutation(self, chainID, residueID, wildtypeAA, mutantAA):
		self.dict["Mutations"].append({
			FieldNames.Chain 		: chainID,
			FieldNames.ResidueID	: residueID,
			FieldNames.WildTypeAA	: wildtypeAA,
			FieldNames.MutantAA		: mutantAA
			})
	
	def addChain(self, chainID):
		self.dict["ExperimentChains"].append(chainID)
	
	def getChains(self):
		return self.dict["ExperimentChains"]
	
	def setMutantIfUnset(self, mutant):
		if not self.dict[FieldNames.Mutant]:
			self.dict[FieldNames.Mutant] = mutant
	
	def addExperimentalScore(self, sourceID, ddG, numMeasurements = 1):
		self.dict["ExperimentalScores"].append({
			FieldNames.SourceID				: sourceID,
			FieldNames.ddG					: ddG,
			FieldNames.NumberOfMeasurements : numMeasurements
			})
	
	def mergeScores(self, maxStdDeviation = 1.0):
		d = self.dict
		
		n = len(d["ExperimentalScores"])
		if n > 1:
			n = float(n)
			sum = 0
			for experimentalResult in d["ExperimentalScores"]:
				if experimentalResult[FieldNames.NumberOfMeasurements] != 1:
					raise Exception("Cannot merge scores when individual scores are from more than one measurement. Need to add logic to do proper weighting.")
				sum += experimentalResult[FieldNames.ddG]
			mean = sum / n
			squaredsum = 0
			for experimentalResult in d["ExperimentalScores"]:
				diff = (experimentalResult[FieldNames.ddG] - mean)
				squaredsum += diff * diff
			variance = squaredsum / n
			d[FieldNames.ScoreVariance] = variance
			stddev = sqrt(variance)
			d["StdDeviation"] = stddev 
			d["WithinStdDeviation"] = stddev <= maxStdDeviation
		else:
			d[FieldNames.ScoreVariance] = 0
			d["WithinStdDeviation"] = True	
	
	def isEligible(self):
		d = self.dict
		if d["WithinStdDeviation"] == None:
			raise Exception("Standard deviation not yet computed.")
		else:
			return d["WithinStdDeviation"]
	
	def __repr__(self):
		d = self.dict
		str = []
		str.append("%s: %s" % (FieldNames.Structure, d[FieldNames.Structure]))
		str.append("%ss: %s" % (FieldNames.Mutant, join(d["Mutants"].keys(), ', ')))
		str.append("%s: %s" % (FieldNames.Source, d[FieldNames.Source]))
		str.append("Chains: %s" % (join([chain for chain in d["ExperimentChains"]], ", ")))
		str.append("Mutations:")
		for mutation in d["Mutations"]:
			str.append("\t%s%d: %s -> %s" % (mutation[FieldNames.Chain], mutation[FieldNames.ResidueID], mutation[FieldNames.WildTypeAA], mutation[FieldNames.MutantAA]))
		str.append("Experimental Scores:")
		for score in d["ExperimentalScores"]:
			n = score[FieldNames.NumberOfMeasurements]
			if n > 1:
				str.append("\t%s\t%0.2f (%d measurements)" % (score[FieldNames.SourceID], score[FieldNames.ddG], score[FieldNames.NumberOfMeasurements]))
			else:
				str.append("\t%s\t%0.2f" % (score[FieldNames.SourceID], score[FieldNames.ddG]))
		return join(str, "\n")
	
	def commit(self, db):
		d = self.dict
		
		if not d[FieldNames.ScoreVariance]:
			self.mergeScores()
		
		if d["Mutants"]:
			for mutant in d["Mutants"].keys():
				MutantStructure = PDBStructure(mutant)
				MutantStructure.commit(db)
			
		#todo self.dbID
		#{FieldNames.Chain	: chainID}
		pass



class Prediction(DBObject):
	
	def __init__(self, ExperimentID, PredictionSet, ToolID, ddG, NumberOfMeasurements = 1):
		self.dict = {
			FieldNames.ExperimentID			: ExperimentID,
			FieldNames.PredictionSet		: PredictionSet,
			FieldNames.ToolID				: ToolID,
			FieldNames.CommandID			: None,
			FieldNames.KeptHETATMLines		: None,
			FieldNames.StrippedPDB			: None,
			FieldNames.InputFiles			: {},
			FieldNames.Description			: {},
			FieldNames.ddG					: ddG,
			FieldNames.NumberOfMeasurements	: NumberOfMeasurements,
		}
		self.dbID = None
		if ExperimentID == None:
			raise Exception("Cannot create the following Prediction - Missing ExperimentID:\n%s\n\t" % self)
		
	def setOptional(self, CommandID = None, KeptHETATMLines = None, StrippedPDB = None, InputFiles = None, Description = None):
		d = self.dict
		if CommandID:
			d[FieldNames.CommandID] = CommandID
		if KeptHETATMLines:
			d[FieldNames.KeptHETATMLines] = KeptHETATMLines
		if StrippedPDB:
			d[FieldNames.StrippedPDB] = StrippedPDB
		if InputFiles:
			d[FieldNames.InputFiles] = InputFiles
		if Description:
			d[FieldNames.Description] = Description			
		
	def getDBID(self):
		if not self.dbID:
			raise Exception("Cannot get the database ID of an uncommitted experiment.")
		else:
			return self.dbID
		
	def commit(self, db):
		#todo
		pass

	def __repr__(self):
		d = self.dict
		str = []
		str.append("%s: %s" % (FieldNames.ExperimentID, d[FieldNames.ExperimentID]))
		str.append("%s: %s" % (FieldNames.PredictionSet, d[FieldNames.PredictionSet]))
		str.append("%s: %d" % (FieldNames.ToolID, d[FieldNames.ToolID]))
		if d[FieldNames.CommandID]:
			str.append("%s: %d" % (FieldNames.CommandID, d[FieldNames.CommandID]))
		if d[FieldNames.KeptHETATMLines] == None:
			str.append("%s: NULL" % (FieldNames.KeptHETATMLines))
		else:
			str.append("%s: %d" % (FieldNames.KeptHETATMLines, d[FieldNames.KeptHETATMLines]))
		n = d[FieldNames.NumberOfMeasurements]
		if n > 1:
			str.append("%s: %0.2f (%d measurements)" % (FieldNames.ddG, d[FieldNames.ddG], n))
		else:
			str.append("%s: %0.2f" % (FieldNames.ddG, d[FieldNames.ddG]))
		
		str.append("%s:" % (FieldNames.InputFiles))
		if d[FieldNames.InputFiles]:
			for k,v in d[FieldNames.InputFiles].iteritems():
				str.append("\t%s" % k)
		else:
			str.append("\tEmpty")
		str.append("%s:" % (FieldNames.Description))
		if d[FieldNames.Description]:
			for k,v in d[FieldNames.Description].iteritems():
				str.append("\t%s: %s" % (k, v))
		else:
			str.append("\tEmpty")
		
		return join(str, "\n")

	def __getitem__(self, key):
		return dict_[key]
	
			
class ddGDatabase(object):
	
	def __init__(self):
		self.connection = MySQLdb.Connection(host = "kortemmelab.ucsf.edu", db = "ddG", user = "kortemmelab", passwd = "r2(#J}(K", port = 3306, unix_socket = "/var/lib/mysql/mysql.sock")
		self.numTries = 32

	def close(self):
		self.connection.close()

	def execQuery(self, sql, parameters = None, cursorClass = MySQLdb.cursors.Cursor, quiet = False):
		"""Execute SQL query"""
		i = 0
		errcode = 0
		caughte = None
		while i < self.numTries:
			i += 1
			try:    
				cursor = self.connection.cursor(cursorClass)
				if parameters:
					errcode = cursor.execute(sql, parameters)
				else:
					errcode = cursor.execute(sql)
				self.connection.commit()
				results = cursor.fetchall()
				cursor.close()
				return results
			except MySQLdb.OperationalError, e:
				errcode = e[0]
				self.connection.ping()
				caughte = e
				continue
			except:                
				traceback.print_exc()
				break
		
		if not quiet:
			sys.stderr.write("\nSQL execution error in query %s at %s:" % (sql, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
			sys.stderr.write("\nErrorcode/Error: %d - '%s'.\n" % (errcode, str(caughte)))
			sys.stderr.flush()
		raise MySQLdb.OperationalError(caughte)

class DatabasePrimer(object):
	'''This class fills in initial values for Tool, AminoAcid, UniProtKB, and UniProtKBMapping. The last will print errors if the corresponding PDB is not in the database.'''
	
	def __init__(self, ddGdb):
		self.ddGdb = ddGdb
		if False:
			self.insertPotapovTools()
			self.insertAminoAcids()
			self.insertUniProtKB()
		
	def insertPotapovTools(self):
		SQL = 'INSERT INTO Tool (Name) VALUES ("CC/PBSA"), ("EGAD"), ("FoldX"), ("Hunter"), ("IMutant2")'
		self.ddGdb.execQuery(SQL)
	
	def insertUniProtKB(self):
		uniprot = os.path.join("rawdata", "uniprotmapping.csv")
		F = open(uniprot)
		lines = F.read().split("\n")[1:]
		F.close()
		UniProtKB = {}
		UniProtKBMapping = {}
		
		for line in lines:
			data = line.split("\t")
			if len(data) == 3:
				PDBID, AC, ID = data 
				UniProtKB[AC] = ID
				UniProtKBMapping[AC] = UniProtKBMapping.get(AC, []) or []
				UniProtKBMapping[AC].append(PDBID)
		
		for AC, ID in sorted(UniProtKB.iteritems()):
			if not self.ddGdb.execQuery("SELECT * FROM UniProtKB WHERE UniProtKB_AC=%s", parameters = (AC,)):
				SQL = 'INSERT INTO UniProtKB (UniProtKB_AC, UniProtKB_ID) VALUES (%s, %s);'
				self.ddGdb.execQuery(SQL, parameters = (AC, ID))
		
		for AC, pdbIDs in sorted(UniProtKBMapping.iteritems()):
			for pdbID in pdbIDs:
				if not self.ddGdb.execQuery("SELECT * FROM UniProtKBMapping WHERE UniProtKB_AC=%s AND PDB_ID=%s", parameters = (AC, pdbID)):
					SQL = 'INSERT INTO UniProtKBMapping (UniProtKB_AC, PDB_ID) VALUES (%s, %s);'
					try:
						self.ddGdb.execQuery(SQL, parameters = (AC, pdbID), quiet = True)
					except:
						pass
			
		
	def insertAminoAcids(self):
		result = self.ddGdb.execQuery("SELECT * FROM AminoAcids")
		print("**")
		print(result)
		print("**")
		aas = [
			["A", "ALA", "Alanine",			"non-polar",	"small"], 
			["C", "CYS", "Cysteine",		"non-polar",	"small"],
			["D", "ASP", "Aspartic acid",	"polar",		"small"],
			["E", "GLU", "Glutamic acid",	"polar",		"large"],
			["F", "PHE", "Phenylalanine",	"non-polar",	"large"],
			["G", "GLY", "Glycine",			"non-polar",	"small"],
			["H", "HIS", "Histidine",		"polar",		"large"],
			["I", "ILE", "Isoleucine",		"non-polar",	"large"],
			["K", "LYS", "Lysine",			"polar",		"large"],
			["L", "LEU", "Leucine",			"non-polar",	"large"],
			["M", "MET", "Methionine",		"non-polar",	"large"],
			["N", "ASN", "Asparagine",		"polar",		"small"],
			["P", "PRO", "Proline",			"non-polar",	"small"],
			["Q", "GLN", "Glutamine",		"polar",		"large"],
			["R", "ARG", "Arginine",		"polar",		"large"],
			["S", "SER", "Serine",			"polar",		"small"],
			["T", "THR", "Threonine",		"polar",		"small"],
			["V", "VAL", "Valine",			"non-polar",	"small"],
			["W", "TRP", "Tryptophan",		"non-polar",	"large"],
			["Y", "TYR", "Tyrosine",		"polar",		"large"]
		]
		for aa in aas:
			SQL = 'INSERT INTO AminoAcids (Code, LongCode, Name, Polarity, Size) VALUES (%s, %s, %s, %s, %s);'
			self.ddGdb.execQuery(SQL, parameters = tuple(aa))

if __name__ == "__main__":
	ddGdb = ddGDatabase()
	DatabasePrimer(ddGdb)
	
