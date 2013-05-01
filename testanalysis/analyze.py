import sys
sys.path.insert(0, "..")
sys.path.insert(0, "../common")
sys.path.insert(0, "../ddglib")
from common import colortext, rosettadb
import analysis

filetype = "pdf"

import ddgdbapi
ddGdb = ddgdbapi.ddGDatabase()
results = ddGdb.execute("SELECT ID, ExperimentID, Status FROM Prediction WHERE PredictionSet=%s", parameters=('AllExperimentsProtocol16',))

single_passed = 0
single_failed = 0
multiple_passed = 0
multiple_failed = 0
for r in results:
	nummutations = len(ddGdb.execute("SELECT * FROM ExperimentMutation WHERE ExperimentID=%s", parameters=(r["ExperimentID"],)))
	if nummutations > 1:
		if r["Status"] == 'done':
			multiple_passed += 1
		elif r["Status"] == 'failed':
			multiple_failed += 1
		else:
			raise Exception("Sd")
	elif nummutations == 1:
		if r["Status"] == 'done':
			single_passed += 1
		elif r["Status"] == 'failed':
			single_failed += 1
		else:
			raise Exception("Sd")
	else:
		raise Exception("Sd")
print("single_passed",single_passed)
print("single_failed",single_failed)
print("multiple_passed",multiple_passed)
print('multiple_failed',multiple_failed)

if False:
	analyzer = analysis.Analyzer("AllExperimentsProtocol16", ddG_score_type = 'kellogg.total')
	analyzer.AddPublishedDDGsToAnalysisTables()
	reporter = analysis.Reporter(analyzer)
	reporter.CreateReport(description = analyzer.description, outfname = 'kellogg.pdf', filetype = filetype)

for score_type in ['noah_8,0A', 'noah_9,0A']:
	for score_method in ['total', 'positional', 'positional_twoscore']:
		analyzer = analysis.Analyzer("AllExperimentsProtocol16", ddG_score_type = '%s.%s' % (score_type, score_method))
		analyzer.AddPublishedDDGsToAnalysisTables()
		reporter = analysis.Reporter(analyzer)
		reporter.CreateReport(description = analyzer.description, outfname = '%s_%s.pdf' % (score_type, score_method), filetype = filetype)
		
#analysis.plot(analysis._R_correlation_coefficient, analysis._createAveragedInputFile, results, "my_plot2.pdf", average_fn = analysis._mean)
		
		
#analysis_objects = analyzer.PlotAll(filetype = "pdf", createFiles = False)
#for o in analysis_objects:
#	print o
#analyzer.plot("Kellogg", analysis.RInterface.correlation_coefficient, "kellogg_R.png", filetype="png")
#analyzer.plot("Kellogg", analysis.RInterface.correlation_coefficient_unfixed, "kellogg_R_unfixed.pdf")
#analyzer.plot("Kellogg", analysis.RInterface._R_mean_absolute_error, "kellogg_R_MAE.pdf")
#analyzer.plot("Potapov", analysis.RInterface.correlation_coefficient, "potapov_R.pdf")
#analyzer.plot("Potapov", analysis.RInterface.correlation_coefficient_unfixed, "potapov_R_unfixed.pdf")
#analyzer.plot("Guerois", analysis.RInterface.correlation_coefficient, "guerois_R.pdf")
#analyzer.plot("Guerois", analysis.RInterface.correlation_coefficient_unfixed, "guerois_R_unfixed.pdf")
#analyzer.plot("ProTherm", analysis.RInterface.correlation_coefficient, "protherm_R.pdf")
#analyzer.plot("ProTherm", analysis.RInterface.correlation_coefficient_unfixed, "protherm_R_unfixed.pdf")


#			"kellogg.txt", "Kellogg")
#for table_name, a_table in sorted(analyzer.analysis_tables.iteritems()):
#	print(a_table)
#	print(table_name)
