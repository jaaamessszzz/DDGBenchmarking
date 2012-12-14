import sys
sys.path.insert(0, "..")
sys.path.insert(0, "../common")
sys.path.insert(0, "../ddglib")
from common import colortext, rosettadb
import pickle
from process import Popen
from common.rosettahelper import readBinaryFile, makeTemp755Directory, writeFile, readFileLines
import zipfile
import shutil
import os
from pdb import PDB, ResidueID2String, checkPDBAgainstMutations, aa1

import ddgdbapi
ddGdb = ddgdbapi.ddGDatabase()
ddGPredictiondb = ddgdbapi.ddGPredictionDataDatabase()

class NoahScore(object):
	
	def __init__(self, total = None, positional = None, positional_twoscore = None):
		self.total = total
		self.positional = positional
		self.positional_twoscore = positional_twoscore
		
	def calculate(self, list_of_files, rosetta_chain, rosetta_resid, radius = 6.0):
		from string import join
		cmd = ([
			'./score_residue.default.linuxgccrelease',
			'--database=/var/binarybuilder/r3.4/rosetta_database',
			'-s']
			+ list_of_files
			+ ['-score:fa_max_dis %0.1f' % radius]
			+ ['-score_residue::residue', '%(rosetta_resid)s%(rosetta_chain)s' % vars()])
		a = join(cmd, " ")
		writeFile("test.sh", a)
		poutput = Popen('.', cmd)
		
		if poutput.errorcode:
			raise Exception("Return code = %s.\n%s" % (str(poutput.errorcode), poutput.stderr))
		else:
			found = {}
			mintotal = None
			minpos = None
			minpostwo = None
			for line in poutput.stdout.split("\n"):
				if line.startswith("apps.score_residue: MINIMUM TOTAL"):
					mintotal = float(line.split()[-1])
				elif line.startswith("apps.score_residue: MINIMUM POSITION (TWOBODY ONLY)"):
					minpostwo = float(line.split()[-1])
				elif line.startswith("apps.score_residue: MINIMUM POSITION"):
					minpos = float(line.split()[-1])
			if not(mintotal != None and minpos != None and minpostwo != None):
				raise Exception("Could not determine all values.")
		self.total = mintotal
		self.positional = minpos
		self.positional_twoscore = minpostwo
	
	def ddg(self, other):
		return NoahScore(
			total = other.total - self.total,
			positional = other.positional - self.positional,
			positional_twoscore = other.positional_twoscore - self.positional_twoscore)
	
	def __repr__(self):
		return("Total score: %(total)f\nPositional score: %(positional)f\nPositional score (two-body) %(positional_twoscore)f" % self.__dict__)

def main():
	results = ddGdb.execute("SELECT ID, ExperimentID, ddG FROM Prediction WHERE PredictionSet='AllExperimentsProtocol16' AND Status='done' AND ScoreVersion <> 0.23")
	if results:
		print("Score versions found which are not 0.23. Need to update table structure.")
		sys.exit(1)
	else:
		# Little hacky way to run two processes
		print( os.getcwd())
		if os.getcwd().endswith("testing_scoring2"):
			results = ddGdb.execute("SELECT ID, ExperimentID, ddG FROM Prediction WHERE PredictionSet='AllExperimentsProtocol16' AND Status='done' AND ScoreVersion = 0.23 AND MOD(ID,2)=1")
		elif os.getcwd().endswith("testing_scoring"):
			results = ddGdb.execute("SELECT ID, ExperimentID, ddG FROM Prediction WHERE PredictionSet='AllExperimentsProtocol16' AND Status='done' AND ScoreVersion = 0.23 AND MOD(ID,2)=0")
		else:
			print("wrong directory")
			sys.exit(0)
	
	import time
	print(len(results))
	#exF = open("existingvalues.txt","w")
	count = 0
	cases_computed = 0
	total_time_in_secs = 0
	
	jobsdone = 0
	jobsleft = 0
				
	for r in results:
		inner_count = 0
		pdbID = ddGdb.execute('SELECT Structure FROM Experiment WHERE ID=%s', parameters=(r['ExperimentID'],))[0]['Structure']
		mutations = ddGdb.execute('SELECT * FROM ExperimentMutation WHERE ExperimentID=%s', parameters=(r['ExperimentID'],))
		extracted_data = False
		
		count += 1
		if len(mutations) == 1:
			timestart = time.time()
			 
			mutation = mutations[0]
			wtaa = mutation['WildTypeAA'] 
			dbchain = mutation['Chain'] 
			mutantaa = mutation['MutantAA']
			
			ddG_dict = pickle.loads(r['ddG'])
			if ddG_dict['version'] != '0.23':
				print(ID, ddG_dict['version'])
			if ddG_dict['version'] == '0.2':
				ddG = ddG_dict['data']['kellogg']['total']['ddG']
			elif ddG_dict['version'] == '0.1':
				ddG = ddG_dict['data']['ddG']
			elif ddG_dict['version'] == '0.21':
				ddG = ddG_dict['data']['kellogg']['total']['ddG']
				ddG_dict['data']['noah_6.0A'] = ddG_dict['data']['noah']
				del ddG_dict['data']['noah']
				ddG_dict['version'] = '0.22'
				pickled_ddG = pickle.dumps(ddG_dict)
				ddGdb.execute('UPDATE Prediction SET ddG=%s WHERE ID=%s', parameters=(pickled_ddG, r['ID'],))
			elif ddG_dict['version'] == '0.22':
				ddG = ddG_dict['data']['kellogg']['total']['ddG']
				try:
					ddG_dict['data']['noah_6,0A'] = ddG_dict['data']['noah_6.0A']
					del ddG_dict['data']['noah_6.0A']
				except:
					pass
				try:
					ddG_dict['data']['noah_7,0A'] = ddG_dict['data']['noah_7.0A']
					del ddG_dict['data']['noah_7.0A']
				except:
					pass
				try:
					ddG_dict['data']['noah_8,0A'] = ddG_dict['data']['noah_8.0A']
					del ddG_dict['data']['noah_8.0A']
				except:
					pass
				try:
					ddG_dict['data']['noah_9,0A'] = ddG_dict['data']['noah_9.0A']
					del ddG_dict['data']['noah_9.0A']
				except:
					pass
				ddG_dict['version'] = '0.23'
				pickled_ddG = pickle.dumps(ddG_dict)
				ddGdb.execute('UPDATE Prediction SET ddG=%s WHERE ID=%s', parameters=(pickled_ddG, r['ID'],))
			elif ddG_dict['version'] == '0.23':
				ddG = ddG_dict['data']['kellogg']['total']['ddG']
			else:
				raise Exception("Booya!")
			
			assert(ddG_dict['version'] == '0.23')
			print(ddG_dict)
			
			all_done = True
			for radius in [7.0, 8.0, 9.0]:
				score_name = ('noah_%0.1fA' % radius).replace(".", ",")
				if not(ddG_dict['data'].get(score_name)):
					all_done = False
				else:
					jobsdone += 1
			if all_done:
				continue
			
			# Extract data
			data = ddGPredictiondb.execute('SELECT Data FROM PredictionData WHERE ID=%s', parameters=(r['ID'],))
			if len(data) == 0:
				colortext.error('No data for id %d' % r['ID'])
				continue
			archivefile = data[0]['Data']
			zipfilename = "%d.zip" % r['ID']
			F = open(zipfilename, "wb")
			F.write(archivefile)
			F.close()
			zipped_content = zipfile.ZipFile(zipfilename, 'r', zipfile.ZIP_DEFLATED)
			tmpdir = None
			repacked_files = []
			mutant_files = []
			try:
				tmpdir = makeTemp755Directory('.')
				highestIndex = -1
				foundResfile = False
				for fname in sorted(zipped_content.namelist()):
					if fname.endswith(".pdb"):
						if fname.startswith("%s/mut_" % r['ID']) or fname.startswith("%s/repacked_" % r['ID']):
							structnum = int(fname[fname.rindex('_')+1:-4]) 
							if fname.startswith("%s/mut_" % r['ID']):
								newfname = 'mutant_%02d' % structnum
							if fname.startswith("%s/repacked_" % r['ID']):
								newfname = 'repacked_%02d' % structnum
							highestIndex = max(highestIndex, structnum)
							
							newfilepath = os.path.join(tmpdir, newfname)
							writeFile(newfilepath, zipped_content.read(fname))
							
							if fname.startswith("%s/mut_" % r['ID']):
								mutant_files.append(newfilepath)
							if fname.startswith("%s/repacked_" % r['ID']):
								repacked_files.append(newfilepath)
						#elif fname.startswith("%s/%s-%s" % (r['ID'],r['ExperimentID'],pdbID)) or fname.startswith("%s/repacked_" % r['ID']):
						#	writeFile(os.path.join(tmpdir, '%s.pdb' % pdbID), zipped_content.read(fname))
					if fname.startswith("%s/%s-%s.resfile" % (r['ID'],r['ExperimentID'],pdbID)):
						foundResfile = True
						lines = zipped_content.read(fname).split("\n")
						assert(len(lines) == 3)
						assert(lines[0] == "NATAA")
						assert(lines[1] == "start")
						resfile_mutation = lines[2].split(" ")
						assert(len(resfile_mutation) == 4)
						rosetta_resid = resfile_mutation[0]
						rosetta_chain = resfile_mutation[1]
						rosetta_mutaa = resfile_mutation[3]
						assert(mutantaa == rosetta_mutaa)
						assert(dbchain == rosetta_chain)
						assert(resfile_mutation[2] == 'PIKAA')
						assert(len(rosetta_mutaa) == 1)
				
				# Make sure the wtaa->mutantaa types match the structures 
				assert(foundResfile)	
				if rosetta_resid.isdigit():
					fullresid = '%s%s%s ' % (rosetta_chain, (4-len(rosetta_resid)) * ' ', rosetta_resid)
				else:
					fullresid = '%s%s%s' % (rosetta_chain, (5-len(rosetta_resid)) * ' ', rosetta_resid)
				assert(PDB(repacked_files[0]).ProperResidueIDToAAMap()[fullresid] == wtaa)
				assert(PDB(mutant_files[0]).ProperResidueIDToAAMap()[fullresid] == mutantaa)
				
				for radius in [7.0, 8.0, 9.0]:
					score_name = ('noah_%0.1fA' % radius).replace(".", ",")
					if ddG_dict['data'].get(score_name):
						continue
					jobsdone += 1

					print("Prediction ID: %d. Calculating radius %0.1f. Calculation #%d of of %d." % (r['ID'], radius, jobsdone, 3*len(results)))
		
					repacked_score = NoahScore()
					repacked_score.calculate(repacked_files, rosetta_chain, rosetta_resid.strip(), radius = radius)
					colortext.message("Repacked")
					print(repacked_score)
					
					mutant_score = NoahScore()
					mutant_score.calculate(mutant_files, rosetta_chain, rosetta_resid.strip(), radius = radius)
					colortext.printf("Mutant", color = 'cyan')
					print(mutant_score)
					
					colortext.printf("ddG", color = 'lightpurple')
					ddg_score = repacked_score.ddg(mutant_score)
					print(ddg_score)
					
					colortext.printf("Liz's ddG", color = 'yellow')
					print("Total score: %f" % ddG)
					
					if ddG_dict['version'] == '0.1':
						ddG_dict['version'] = '0.21'
						ddG_dict['data'] = {
							'kellogg' : {
								'total' : ddG_dict['data'],
							},
							'noah': {
								'total' : {'ddG' : ddg_score.total},
								'positional' : {'ddG' : ddg_score.positional},
								'positional_twoscore' : {'ddG' : ddg_score.positional_twoscore},
							},
						}
					elif ddG_dict['version'] == '0.2':
						ddG_dict['version'] = '0.21'
						ddG_dict['data']['noah']['total']['ddG'] = ddg_score.total 
						ddG_dict['data']['noah']['positional']['ddG'] = ddg_score.positional 
						ddG_dict['data']['noah']['positional_twoscore']['ddG'] = ddg_score.positional_twoscore 
					elif ddG_dict['version'] == '0.22':
						ddG_dict['data'][score_name] = {'total' : {}, 'positional' : {}, 'positional_twoscore' : {}}
						ddG_dict['data'][score_name]['total']['ddG'] = ddg_score.total 
						ddG_dict['data'][score_name]['positional']['ddG'] = ddg_score.positional 
						ddG_dict['data'][score_name]['positional_twoscore']['ddG'] = ddg_score.positional_twoscore 
					elif ddG_dict['version'] == '0.23':
						ddG_dict['data'][score_name] = {'total' : {}, 'positional' : {}, 'positional_twoscore' : {}}
						ddG_dict['data'][score_name]['total']['ddG'] = ddg_score.total 
						ddG_dict['data'][score_name]['positional']['ddG'] = ddg_score.positional 
						ddG_dict['data'][score_name]['positional_twoscore']['ddG'] = ddg_score.positional_twoscore 
						
					pickled_ddG = pickle.dumps(ddG_dict)
					ddGdb.execute('UPDATE Prediction SET ddG=%s WHERE ID=%s', parameters=(pickled_ddG, r['ID'],))
				shutil.rmtree(tmpdir)
				os.remove(zipfilename)
					
			except Exception, e:
				print("Exception!", str(e))
				import traceback
				print(traceback.format_exc())
				if tmpdir:
					shutil.rmtree(tmpdir)

			cases_computed += 1
			timetaken_in_secs = time.time() -timestart
			total_time_in_secs += timetaken_in_secs
			number_of_cases_left = len(results) - count
			average_time_taken = float(total_time_in_secs)/float(cases_computed)
			estimate_remaining_time = number_of_cases_left * average_time_taken
			
			colortext.message("Time taken for this case: %0.2fs." % timetaken_in_secs)
			colortext.message("Average time taken per case: %0.2fs." % average_time_taken)
			colortext.message("Estimated time remaining: %dh%dm%ds." % (int(estimate_remaining_time/3600), int((estimate_remaining_time/60) % 60), estimate_remaining_time % 60))
			print("\n")
	print('jobsdone', jobsdone)
	print('jobsleft', jobsleft)
	
	#exF.close()
	
main()
	
