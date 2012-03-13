import sys
import os
import re
import pickle
import subprocess
from string import join
import time
sys.path.insert(0, "../common")
import colortext
import ddgproject
from tempfile import mkstemp

ddGdb = ddgproject.ddGDatabase()
dbfields = ddgproject.FieldNames()

def _mean(points, numpoints):
	'''Points is expected to be a list (or iterable), numpoints should be an integer.'''
	return sum(points) / numpoints 

def _createMAEFile(results, outfname, average_fn = _mean):
	output = []
	output.append("X,Error")
	X = 1
	for r in results:
		predicted_ddG = pickle.loads(r["ddG"])["data"]["ddG"]
		
		expID = r["ExperimentID"]
		expScores = ddGdb.callproc('GetScores', parameters = (expID,))
		
		n = 0
		scores = []
		for e in expScores:
			c = e["NumberOfMeasurements"]
			n += c
			scores.append(e["ddG"] * c)
		point = abs(average_fn(scores, n) - predicted_ddG)
		output.append("%s,%s"% (X,point))
		X += 1
		
	F, fname = mkstemp(dir = ".")
	F = os.fdopen(F, "w")
	F.write(join(output, "\n"))
	F.write("\n")
	F.close()
	return fname

def _createAveragedInputFile(results, outfname, average_fn = _mean):
	output = []
	output.append("Experimental,Predicted")
	todo = 0
	for r in results:
		predicted_ddG = pickle.loads(r["ddG"])["data"]["ddG"] + todo
		todo += 0.23
		expID = r["ExperimentID"]
		expScores = ddGdb.callproc('GetScores', parameters = (expID,))
		
		n = 0
		scores = []
		for e in expScores:
			c = e["NumberOfMeasurements"]
			n += c
			scores.append(e["ddG"] * c)
		eavg = average_fn(scores, n)
		
		output.append("%s,%s"% (eavg, predicted_ddG))

	F, fname = mkstemp(dir = ".")
	F = os.fdopen(F, "w")
	F.write(join(output, "\n"))
	print(join(output, "\n"))
	F.write("\n")
	F.close()
	return fname

def _runRScript(RScript):
	F, rscriptname = mkstemp(dir = ".")
	F = os.fdopen(F, "w")
	F.write(RScript)
	F.close()
	p = subprocess.Popen(["R","CMD", "BATCH", rscriptname]) 
	while True:
		time.sleep(0.3)
		errcode = p.poll()
		if errcode != None:
			break
	rout = "%s.Rout" % rscriptname
	os.remove(rscriptname)
	if errcode != 0:
		if os.path.exists(rout):
			F = open(rout, "r")
			colortext.warning(F.read())
			F.close()
			os.remove(rout)
		raise colortext.Exception("The R script failed with error code %d." % errcode)
	os.remove(rout)
	
def _R_mean_unsigned_error(inputfname, outfname):
	RScript='''
pdf('%(outfname)s')
par(mar=c(5, 5, 1, 1))
a <- read.csv('%(inputfname)s', header=T)
head(a)
plot(a$Error, pch=19,xlab='Experiment #',ylab='MAE',cex.lab=1.5, xlim=c(1,max(a$X)+1), ylim=c(0,max(a$Error)+1), xaxt='n')
axis(side=1, at=0:23, cex.axis=.5)
?axis
dev.off()''' % vars()
	# todo: The graph size is not being set correctly i.e. the max function call is wrong above
	_runRScript(RScript)


def _R_correlation_coefficient(inputfname, outfname):
	RScript='''
pdf('%(outfname)s')
par(mar=c(5, 5, 1, 1))
a <- read.csv('%(inputfname)s', header=T)
head(a)
reg1 <- lm(a$Experimental~a$Predicted)
plot(a$Experimental, a$Predicted, pch=19,xlab='Experimental ddG',ylab='Computational ddG',cex.lab=1.5, xlim=c(0,max(a$Predicted)+3), ylim=c(0,max(a$Experimental)+3))
abline(reg1)
dev.off()''' % vars()
	# todo: The graph size is not being set correctly i.e. the max function call is wrong above
	_runRScript(RScript)


def plot(RFunction, filecreator, results, outfname, average_fn = _mean):
	'''Results is expect to be a list of dicts each of which has the keys ExperimentID and ddG.'''
	if not results:
		raise Exception("The results set is empty.")
	elif len(results) == 1:
		raise Exception("The results set only has one data point. At least two points are required.")
	else:
		inputfname = filecreator(results, outfname, average_fn)
		try:
			RFunction(inputfname, outfname)
		except Exception, e:
			import traceback
			colortext.error(traceback.format_exc())
			os.remove(inputfname)
			raise Exception(e)
		os.remove(inputfname)
			

#results = ddGdb.execute('SELECT ExperimentID, ddG FROM Prediction WHERE PredictionSet="testrun" AND ddG IS NOT NULL')
#t1 = {"ExperimentID" : (results[0]["ExperimentID"] + 1), "ddG" : (results[0]["ddG"])}
#t2 = {"ExperimentID" : (results[0]["ExperimentID"] + 2), "ddG" : (results[0]["ddG"])}
#results = (results[0], t1 ,t2) 
#plot(_R_mean_unsigned_error, _createMAEFile, results, "my_plot1.pdf", average_fn = _mean)
#plot(_R_correlation_coefficient, _createAveragedInputFile, results, "my_plot2.pdf", average_fn = _mean)

# par is used to specify graphical parameters. mar sets the margin sizes - c(bottom, left, top right)
# read.table is useful for reading delimited files. e.g. read.table("c:/mydata.csv", header=TRUE, sep=",", row.names="id")
# plot creates a plot
#	pch=19 means use a diamond character (web-search for 'R plot symbols' pch)
#	a$V1, a$V2 choose the first and second columns of data table a
#	xlab and ylab set the labels
#	cex	indicates the amount by which plotting text and symbols should be scaled relative to the default. 1=default, 1.5 is 50% larger, 0.5 is 50% smaller, etc.
#	cex.lab	defines the magnification of x and y labels relative to cex
#	abline(k,l) adds a straight line where k and l are the intercept and slope values respectively








class kobject(object):
	
	def isOfClass(self, c):
		return self.__class__ == c

	def getBaseClassesForObject(self):
		return self.getBaseClasses()

	def getClassName(self):
		return self.__class__.__name__

	@staticmethod
	def _getBaseClasses(c, b):
		superclasses = [sc for sc in c.__bases__ if sc not in b]
		b += superclasses
		for sc in superclasses:
			kobject._getBaseClasses(sc, b)
		return b
	
	@classmethod
	def getBaseClasses(cls):
		return kobject._getBaseClasses(cls, [])
	
class Filter(kobject):
	
	_ResultSet = None
	
	def __init__(self):
		self.paramsoffset = 0
		#raise Exception("Calling an abstract class.")
		
	def apply(self, db):
		self.conditions = []
		self.parameters = []
		
		self._apply()
		
		conditions = self.conditions
		parameters = self.parameters
		assert(len(conditions) == len(parameters) - self.paramsoffset)
		if conditions:
			result_set = self._ResultSet(db, SQL = "WHERE %s" % join(conditions, " AND "), parameters = tuple(parameters))
		else:
			result_set = self._ResultSet(db)

		for postfilter in self.post_SQL_filters:
			# post_SQL_filters functions take a database connection and self._ResultSet object and return a self._ResultSet object 
			result_set = postfilter(db, result_set)
		
		return result_set

	def __or__(self, disjunct):
		if disjunct.isOfClass(UnionFilter):
			return UnionFilter([self] + disjunct.filters)
		elif Filter in disjunct.getBaseClasses():
			return UnionFilter([self, disjunct])
		else:
			raise Exception("Trying to combine a filter with an object of class '%s'." % disjunct.getClassName())

class UnionFilter(kobject):
	
	def __init__(self, filters = []):
		self.filters = filters
	
	def add(self, filter):
		if Filter in filter.getBaseClasses():
			self.filters.append(filter)
		else:
			raise Exception("Trying to combine a filter with an object of class '%s'." % disjunct.getClassName())
	
	def getFilters(self):
		return self.filters
	
	def __or__(self, disjunct):
		if disjunct.isOfClass(UnionFilter):
			return UnionFilter(self.filters + disjunct.filters)
		elif Filter in disjunct.getBaseClasses():
			return UnionFilter(self.filters + [disjunct])
		else:
			raise Exception("Trying to combine a filter with an object of class '%s'." % disjunct.getClassName())
			
class ResultSet(kobject):
	dbname = None
	primary_key = "ID"
	stored_procedures = [] # e.g. ["GetScores"]
	allowed_filters = []
	
	stored_procedure_regex = re.compile("CALL (\w+)")
	
	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = []):
		'''e.g. ResultSet("CALL GetScores", parameters = 16734)
				ResultSet("WHERE ...")'''
		
		self._log = []
		self.db = db
		if parameters:
			if type(parameters) != tuple:
				parameters = (parameters,)
			if not SQL:
				raise Exception("Parameters %s were specified for an empty SQL query." % parameters)
		
		results = []
		if SQL or (not AdditionalIDs):
			# We run an SQL query if the SQL parameter has been specified or if AdditionalIDs is empty.
			# If AdditionalIDs is not empty and the SQL is blank then we do NOT run the SQL query and
			# instead use AdditionalIDs as the record keys.  
			sp = ResultSet.stored_procedure_regex.match(SQL)
			if sp:
				if not sp.groups(1) in self.__class__.stored_procedures:
					raise Exception("Cannot call stored procedure %s to create a %s" % (sp.groups(1)[0], self.__class__.__name__))
				else:
					results = db.callproc(SQL, parameters)
					self.log("ResultSet object initialized with stored procedure call %s%s." % (SQL, parameters))
			else:
				if not self.__class__.dbname:
					raise Exception("There is not database table associated with the class %s." % ( self.__class__.__name__))
				SQL = "SELECT %s FROM %s %s" % (self.__class__.primary_key, self.__class__.dbname, SQL)
				results = self.db.execute(SQL, parameters)
				self.log("ResultSet object initialized with SQL query '%s' %% %s." % (SQL, parameters or ""))
		
		AdditionalIDs = set(AdditionalIDs)
		if AdditionalIDs:
			AdditionalIDs = set(AdditionalIDs)
			allIDs = self.db.execute("SELECT %s FROM %s" % (self.__class__.primary_key, self.__class__.dbname), parameters)
			pIDs = set([r[self.__class__.primary_key] for r in allIDs]).intersection(AdditionalIDs)
			if not (len(pIDs) == len(AdditionalIDs)):
				raise Exception("Records associated with the following IDs could not be found: %s." % join(AdditionalIDs.difference(pIDs), ","))
		
		if SQL:
			if not results:
				raise Exception("No results were returned from '%s' %% %s." % (SQL, parameters or ""))
			else:
				if not self.__class__.primary_key in results[0]:
					raise Exception("The resulting set from '%s'(%s) must including the primary key field %s of %s." % (SQL, parameters, self.__class__.primary_key, self.__class__.dbname))
				
		self.IDs = set([r[self.__class__.primary_key] for r in results]).union(AdditionalIDs)
		self.initialresults = None
		self.filterchain = []
		self.log("Initial record count: %d" % len(self.IDs))
	
	def log(self, msg):
		self._log.append(msg)
				
	def __repr__(self):
		return join(self._log, "\n")
	
	def getInitialResults(self):
		if not self.initialresults:
			SQL = "SELECT * FROM %s" % self.__class__.dbname
			results = self.db.execute(SQL)
			self.initialresults = [r for r in results if r[self.__class__.primary_key] in self.IDs]
		return self.initialresults

	def addFilter(self, filter, tag = None):
		if not (Filter in filter.getBaseClassesForObject() or filter.isOfClass(UnionFilter)):
			raise Exception("Trying to add a filter which does not derive from the Filter class.")
		
		allowed = False 
		for f in self.__class__.allowed_filters:
			if filter.isOfClass(f):
				allowed = True
				break
		if not allowed:
			if filter.isOfClass(UnionFilter):
				for f in filter.getFilters():
					if not (f.__class__ in self.__class__.allowed_filters):
						raise Exception("Trying to add a %s filter which is not supported by this class (%s)." % (f.getClassName(), self.getClassName()))
			else:
				raise Exception("Trying to add a %s filter which is not supported by this class (%s)." % (filter.getClassName(), self.getClassName()))
		
		self.filterchain.append((filter, tag))
	
	def applyUnionFilter(self, pks, ufilter, tag):
		pkset = set()
		for filter in ufilter.getFilters():
			self.log("  Applying %s:" % (tag or filter.getClassName()))
			pkset = pkset.union(self.applyFilter(pks, filter, tag))
		return pkset
	
	def applyFilter(self, pks, filter, tag):
		self.log("Applying %s:" % (tag or filter.getClassName()))
		raise Exception("This function needs to be implemented for this class (%s)." % self.getClassName())
		
	def getFilteredResults(self):
		'''Applies the filters, returns a new dict'''
		pks = self.IDs
		self.log("Primary keys before filtering:")
		self.log(str(pks))
		for taggedFilter in self.filterchain:
			filter = taggedFilter[0] 
			tag = taggedFilter[1]
			if filter.isOfClass(UnionFilter):
				self.log("Applying %s:"% (tag or "Union Filter"))
				pks = self.applyUnionFilter(pks, filter, tag)
				self.log("Primary keys after filtering:")
				self.log(str(pks))
			elif Filter in filter.getBaseClasses():
				self.log("Applying %s:" % (tag or filter.getClassName()))
				pks = self.applyFilter(pks, filter, tag)
				self.log("Primary keys after filtering:")
				self.log(str(pks))
			else:
				raise Exception("BLARG!")
		self.log("Filtered record count: %d" % len(pks))
		
		SQL = "SELECT * FROM %s" % self.__class__.dbname
		results = self.db.execute(SQL)
		return [r for r in results if r[self.__class__.primary_key] in pks]
	
	def intersection(self, IDs):
		return self.IDs.intersection(IDs)

class StructureResultSet(ResultSet):
	dbname = dbfields.Structure
	primary_key = dbfields.PDB_ID
	stored_procedures = []
	allowed_filters = []

	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = []):
		super(StructureResultSet, self).__init__(db, SQL, parameters, AdditionalIDs)
	
	def applyFilter(self, pks, filter, tag):
		if filter.isOfClass(StructureFilter):
			structures = filter.apply(self.db)
			return pks.intersection(structures.IDs)

class StructureFilter(Filter):
	'''Note: minresolution and maxresolution are soft bounds i.e. if they are set to 1 and 2.5 respectively then we find structures with resolution r where 1 <= r <= 2.5.
			 If you want to specify 1 <= r < 2.5 instead, passing 2.49999 as the max resolution should usually work.'''
	
	_ResultSet = StructureResultSet
	
	XRay = "X-RAY DIFFRACTION"
	NMR = "SOLUTION NMR"
	_AllowedTechniques = [XRay, NMR]
	
	# Constructors
	
	def __init__(self):
		super(StructureFilter, self).__init__()
		self.nullresolution = None
		self.minresolution = None
		self.maxresolution = None
		self.bfactors_min = None
		self.bfactors_max = None
		self.techniques = None
		self.post_SQL_filters = []
		
	@staticmethod
	def TotalBFactors(min, max):
		sf = StructureFilter()
		sf.setTotalBFactors(min, max)
		return sf

	@staticmethod
	def WithNullResolution(allowNull = True):
		sf = StructureFilter()
		sf.nullresolution = allowNull
		return sf

	@staticmethod
	def Resolution(min, max):
		sf = StructureFilter()
		sf.setResolution(min, max)
		return sf
	
	@staticmethod
	def Techniques(techniques):
		sf = StructureFilter()
		sf.setTechniques(techniques)
		return sf
	
	# Property setters
	
	def setAllowNullResolution(self, allowNull):
		# This function is added to keep in line with the reset of the API
		sf.nullresolution = allowNull
	
	def setTechniques(self, techniques):
		if type(techniques) != list:
			techniques = [techniques]
		for t in techniques:
			if not t in self._AllowedTechniques:
				raise Exception("Technique %s is not recognized." % t) 
		self.techniques = techniques
	
	def setResolution(self, min, max):
		self.minresolution = min
		self.maxresolution = max
	
	def setTotalBFactors(self, min, max):
		self.bfactors_min = min
		self.bfactors_max = max

	# _apply and filter functions
		
	def _apply(self):
		self.post_SQL_filters = []
		if self.nullresolution != None:
			if self.nullresolution:
				self.conditions.append("Resolution IS NULL")
			else:
				self.conditions.append("Resolution IS NOT NULL")
			self.paramsoffset -=  1
		if self.minresolution != None:
			self.conditions.append("Resolution >= %s")
			self.parameters.append(self.minresolution)
		if self.maxresolution != None:
			self.conditions.append("Resolution <= %s")
			self.parameters.append(self.maxresolution)
		if self.techniques:
			tstrs = []
			for t in self.techniques:
				tstrs.append('Techniques LIKE "%%%s%%"' % t)
			self.conditions.append(join(tstrs, " OR "))
			self.paramsoffset -=  1
		
		if self.bfactors_min or self.bfactors_max:
			self.post_SQL_filters.append(self._checkTotalBFactorRange) 
			
	def _checkTotalBFactorRange(self, db, result_set):
		m_min = self.bfactors_min
		m_max = self.bfactors_max
		pkname = self._ResultSet.primary_key
		
		new_IDs = []
		for record in result_set.getInitialResults():
			bfactors =  pickle.loads(record["BFactors"])
			average = bfactors["Total"][0]
			passed = True
			if m_min and average < m_min: 
				passed = False
			if m_max and average > m_max: 
				passed = False
			if passed:
				new_IDs.append(record[pkname])
		return self._ResultSet(db, AdditionalIDs = new_IDs)
		
StructureResultSet.allowed_filters = [StructureFilter]	
		
class PredictionResultSet(ResultSet):
	dbname = dbfields.Prediction
	primary_key = dbfields.ID
	stored_procedures = []
	allowed_filters = [StructureFilter]
	
	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = []):
		super(PredictionResultSet, self).__init__(db, SQL, parameters, AdditionalIDs)
		
		self.structure_map = {} 
		for id in self.IDs:
			results = db.execute("SELECT Structure.PDB_ID FROM Prediction INNER JOIN Experiment ON ExperimentID=Experiment.ID INNER JOIN Structure on Experiment.Structure=Structure.PDB_ID WHERE Prediction.ID=%s", parameters=(id,), cursorClass=ddgproject.StdCursor)
			pdbID = results[0][0]
			self.structure_map[pdbID] = self.structure_map.get(pdbID) or []
			self.structure_map[pdbID].append(id)
		
	
	def applyFilter(self, pks, filter, tag):
		if filter.isOfClass(StructureFilter):
			structures = filter.apply(self.db)
			foundIDs = []
			for s in structures.IDs:
				if self.structure_map.get(s):
					foundIDs.extend(self.structure_map[s])
			return pks.intersection(foundIDs)

class ExperimentResultSet(ResultSet):
	dbname = dbfields.Experiment
	primary_key = dbfields.ID
	stored_procedures = ["GetScores"]
	allowed_filters = [StructureFilter]
	
	
