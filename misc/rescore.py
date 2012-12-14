import sys
sys.path.insert(0, "..")
sys.path.insert(0, "../common")
sys.path.insert(0, "../ddglib")
import zipfile
import shutil
import os
import time
import pickle
import subprocess
from process import Popen
from common import colortext, rosettadb
from common.rosettahelper import readBinaryFile, makeTemp755Directory, writeFile, readFileLines
from pdb import PDB, ResidueID2String, checkPDBAgainstMutations, aa1
	
import ddgdbapi
ddGdb = ddgdbapi.ddGDatabase()
ddGPredictiondb = ddgdbapi.ddGPredictionDataDatabase()

current_score_revision = '0.23'
class WrongScoreRevisionException(Exception): pass

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

def delete_scores(results, score_type):
	'''e.g. delete_scores(results, 'noah_6,0A')'''
	for r in results:
		if len(mutations) == 1:
			ddG_dict = pickle.loads(r['ddG'])
			if ddG_dict['version'] != current_score_revision:
				raise WrongScoreRevisionException("Expected score revision %s. Found revision %s instead." % (current_score_revision, ddG_dict['version']))
			if ddG_dict['data'].get(score_type):
				del ddG_dict['data'][score_type]
				pickled_ddG = pickle.dumps(ddG_dict)
				ddGdb.execute('UPDATE Prediction SET ddG=%s WHERE ID=%s', parameters=(pickled_ddG, r['ID'],))

def rename_score_type(results, old_score_type, new_score_type):
	'''e.g. rename_score_type(results, 'noah_6,0A', 'noah_6A')'''
	for r in results:
		if len(mutations) == 1:
			ddG_dict = pickle.loads(r['ddG'])
			if ddG_dict['version'] != current_score_revision:
				raise WrongScoreRevisionException("Expected score revision %s. Found revision %s instead." % (current_score_revision, ddG_dict['version']))
			if ddG_dict['data'].get(old_score_type):
				if ddG_dict['data'].get(new_score_type):
					raise WrongScoreRevisionException("Found both the old score type %s and the new score type %s." % (old_score_type, new_score_type))
				ddG_dict['data'][new_score_type] = ddG_dict['data'][old_score_type]
				del ddG_dict['data'][old_score_type]
				pickled_ddG = pickle.dumps(ddG_dict)
				ddGdb.execute('UPDATE Prediction SET ddG=%s WHERE ID=%s', parameters=(pickled_ddG, r['ID'],))
	
def update_records(results):
	print("Updating scores to revision %s." % current_score_revision)
	for r in results:
		inner_count = 0
		extracted_data = False
		
		count += 1
		if len(mutations) == 1:
			ddG_dict = pickle.loads(r['ddG'])
			if ddG_dict['version'] != current_score_revision:
				# update code goes here
				#ddG_dict = pickle.loads(r['ddG'])
				#do something
				#pickled_ddG = pickle.dumps(ddG_dict)
				#ddGdb.execute('UPDATE Prediction SET ddG=%s WHERE ID=%s', parameters=(pickled_ddG, r['ID'],))
				raise WrongScoreRevisionException("Expected score revision %s. Found revision %s instead." % (current_score_revision, ddG_dict['version']))

from optparse import OptionParser

def get_number_of_processors():
	'''Only runs on Linux and I'm unsure how general this is. Works for our webserver.'''
	
	poutput = subprocess.Popen(['more', '/proc/cpuinfo'], stdout=subprocess.PIPE)
	poutput = subprocess.Popen(['grep', 'processor'], stdin=poutput.stdout, stdout=subprocess.PIPE)
	if poutput.returncode:
		raise Exception("Return code = %s.\n%s" % (str(poutput.errorcode), poutput.stderr))
	else:
		return len(poutput.stdout.readlines())

class Timer(object):
	'''Simple class for non-nested timers.'''
	
	max_name_length = 35
	
	def __init__(self):
		self.stages = []
		
	def add(self, stage):
		stage = stage[:Timer.max_name_length]
		if self.stages:
			laststage = self.stages[-1]
			laststage['time_taken'] = time.time() - laststage['timer_start']
			laststage['timer_start'] = None
		self.stages.append({'name' : stage, 'time_taken' : 0, 'timer_start' : time.time()})
	
	def stop(self):
		if self.stages:
			laststage = self.stages[-1]
			laststage['time_taken'] = time.time() - laststage['timer_start']
			laststage['timer_start'] = None
	
	def sum(self):
		return sum([s['time_taken'] for s in self.stages])
	
	def __repr__(self):
		s = []
		for stage in self.stages:
			s.append('%s %fs' % (stage['name'].ljust(Timer.max_name_length + 2), stage['time_taken']))
		return "\n".join(s)
			
def main():
	max_processors = get_number_of_processors() 
	
	rescore_process_file = "/tmp/klab_rescore.txt"
	parser = OptionParser()
	parser.add_option("-n", "--numprocesses", default=1, type='int', dest="num_processes", help="The number of processes used for the rescoring. The cases are split according to this number.", metavar="NUM_PROCESSES")
	parser.add_option("-p", "--process", default=1, type='int', dest="process", help="The ID of this process. This should be an integer between 1 and the number of processes used for the rescoring.", metavar="PROCESS_ID")
	parser.add_option("-d", "--delete",  action="store_true", dest="delete", help="Delete the process tracking file %s." % rescore_process_file)
	(options, args) = parser.parse_args()
	
	if options.delete and os.path.exists(rescore_process_file):
		print("Removing %s." % rescore_process_file)
		os.remove(rescore_process_file)
	
	num_processes = options.num_processes
	process_id = options.process
	if num_processes < 1:
		raise colortext.Exception("At least 1 processor must be used.")
	if num_processes > max_processors:
		raise colortext.Exception("Only %d processors/cores were detected. Cannot run with %d processes." % (max_processors, num_processes))
	if num_processes > (max_processors * 0.75):
		colortext.warning("Warning: Using %d processors/cores out of %d which is %0.2f%% of the total available." % (num_processes, max_processors, (100.0*float(num_processes)/float(max_processors))))
	if not(1 <= process_id <= min(max_processors, num_processes)):
		raise colortext.Exception("The process ID %d must be between 1 and the number of processes, %d." % (process_id, num_processes))
	
	if os.path.exists(rescore_process_file):
		lines = readFileLines(rescore_process_file)
		idx = lines[0].find("numprocesses")
		if idx == -1:
			raise Exception("Badly formatted %s." % rescore_process_file)
		existing_num_processes = int(lines[0][idx+len("numprocesses"):])
		if existing_num_processes != num_processes:
			raise colortext.Exception("You specified the number of processes to be %d but %s already specifies it as %d." % (num_processes, rescore_process_file, existing_num_processes))
		for line in [line for line in lines[1:] if line.strip()]:
			idx = line.find("process")
			if idx == -1:
				raise colortext.Exception("Badly formatted %s. Line is '%s'." % (rescore_process_file, line))
			existing_process = int(line[idx+len('process'):])
			if process_id == existing_process:
				raise colortext.Exception("Process %d is already logged as running. Check if this is so and edit %s." % (process_id, rescore_process_file))
		F = open(rescore_process_file, 'a')
		F.write("process %d\n" % process_id)
		F.close()
	else:
		F = open(rescore_process_file, 'w')
		F.write("numprocesses %d\n" % num_processes)
		F.write("process %d\n" % process_id)
		F.close()
	
	output_dir = 'rescoring/%d' % process_id
	if not(os.path.exists(output_dir)):
		os.makedirs(output_dir)
	abs_output_dir = os.path.abspath(os.path.join(os.getcwd(), output_dir))
	print("Running process in %s.\n" % abs_output_dir)
	
	results = ddGdb.execute("SELECT ID, ExperimentID, ddG FROM Prediction WHERE PredictionSet='AllExperimentsProtocol16' AND Status='done' AND ScoreVersion <> %s", parameters=(float(current_score_revision),))
	if results:
		raise WrongScoreRevisionException("Score versions found which are not %s. Need to update table structure." % current_score_revision)
	else:
		# Hacky way to run multiple processes
		results = ddGdb.execute("SELECT ID, ExperimentID, ddG FROM Prediction WHERE PredictionSet='AllExperimentsProtocol16' AND Status='done' AND ScoreVersion=%s AND MOD(ID,%s)=%s", parameters=(float(current_score_revision),num_processes,process_id-1))
	
	radii = [7.0, 8.0, 9.0]
	
	count = 0
	cases_computed = 0
	total_time_in_secs = 0
	
	number_of_cases_left = len(results) * len(radii)
	
	colortext.printf("Rescoring %d predictions over %d radii...\n" % (len(results), len(radii)), 'lightgreen')
	for r in results:
		t = Timer()
		t.add('Preamble')
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
			kellogg_ddG = ddG_dict['data']['kellogg']['total']['ddG']
			assert(ddG_dict['version'] == current_score_revision)
			
			all_done = True
			for radius in radii:
				score_name = ('noah_%0.1fA' % radius).replace(".", ",")
				if not(ddG_dict['data'].get(score_name)):
					all_done = False
				else:
					cases_computed += 1
					number_of_cases_left -= 1
			if all_done:
				print('Prediction %d: done.' % r["ID"])
				continue
			
			# Extract data
			t.add('Grab data')
			data = ddGPredictiondb.execute('SELECT Data FROM PredictionData WHERE ID=%s', parameters=(r['ID'],))
			if len(data) == 0:
				colortext.error('No data for id %d' % r['ID'])
				continue
			archivefile = data[0]['Data']
			zipfilename = os.path.join(output_dir, "%d.zip" % r['ID'])
			F = open(zipfilename, "wb")
			F.write(archivefile)
			F.close()
			
			t.add('Extract data')
			zipped_content = zipfile.ZipFile(zipfilename, 'r', zipfile.ZIP_DEFLATED)
			tmpdir = None
			repacked_files = []
			mutant_files = []
			try:
				tmpdir = makeTemp755Directory(output_dir)
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
				
				for radius in radii:
					score_name = ('noah_%0.1fA' % radius).replace(".", ",")
					
					if ddG_dict['data'].get(score_name):
						print('Radius %0.1f: done.' % radius)
						continue
					cases_computed += 1
					number_of_cases_left -= 1
					
					t.add('Radius %0.3f: repacked' % radius)
					colortext.printf("Prediction ID: %d. Calculating radius %0.1f. Calculation #%d of %d." % (r['ID'], radius, cases_computed, len(results) * len(radii)), 'orange')
		
					repacked_score = NoahScore()
					repacked_score.calculate(repacked_files, rosetta_chain, rosetta_resid.strip(), radius = radius)
					colortext.message("Repacked")
					print(repacked_score)
					
					t.add('Radius %0.3f: mutant' % radius)
					mutant_score = NoahScore()
					mutant_score.calculate(mutant_files, rosetta_chain, rosetta_resid.strip(), radius = radius)
					colortext.printf("Mutant", color = 'cyan')
					print(mutant_score)
					
					t.add('Radius %0.3f: postamble' % radius)
					colortext.printf("ddG", color = 'lightpurple')
					ddg_score = repacked_score.ddg(mutant_score)
					print(ddg_score)
					
					colortext.printf("Liz's ddG", color = 'yellow')
					print("Total score: %0.3f" % kellogg_ddG)
					
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
				t.add('Cleanup')
				shutil.rmtree(tmpdir)
				os.remove(zipfilename)
					
			except Exception, e:
				print("Exception!", str(e))
				import traceback
				print(traceback.format_exc())
				if tmpdir:
					shutil.rmtree(tmpdir)

			total_time_in_secs += t.sum()
			average_time_taken = float(total_time_in_secs)/float(cases_computed)
			estimate_remaining_time = number_of_cases_left * average_time_taken
			
			t.stop()
			colortext.printf("**Profile**", 'orange')
			print(t)
			colortext.message("Time taken for this case: %0.2fs." % t.sum())
			colortext.message("Average time taken per case: %0.2fs." % average_time_taken)
			colortext.message("Estimated time remaining: %dh%dm%ds." % (int(estimate_remaining_time/3600), int((estimate_remaining_time/60) % 60), estimate_remaining_time % 60))
			print("\n")
			
	#exF.close()
	colortext.printf("\nDone.", 'lightgreen')
	
main()
	
