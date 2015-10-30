import sys
sys.path.insert(0, "..")
sys.path.insert(0, "../common")
sys.path.insert(0, "../ddglib")
sys.path.insert(0, "../..")
import analysis


filetype = "pdf"

prediction_set_for_analysis = 'RosettaCon2013_P16_talaris2013'
prediction_set_for_analysis = 'RosettaCon2013_P16_score12prime'

prediction_set_for_analysis = 'RosCon2013_P16_score12prime'
prediction_set_for_analysis = 'RosCon2013_P16_talaris2013'
prediction_set_for_analysis = 'RosCon2013_P16_score12prime'
prediction_set_for_analysis = 'RosCon2013_P16_talaris2013sc'
prediction_set_for_analysis = 'Protocol_16_r57471'

import ddgdbapi
if __name__ == '__main__':
    ddGdb = ddgdbapi.ddGDatabase()

results = ddGdb.execute("SELECT ID, ExperimentID, Status FROM Prediction WHERE PredictionSet=%s AND Status IN ('done','failed')", parameters=(prediction_set_for_analysis,))

if False:
    FailedJobs = {}
    results = ddGdb.execute("SELECT ID, ExperimentID, UserDataSetExperimentID FROM Prediction WHERE PredictionSet='RosCon2013_P16_score12prime' AND Status IN ('failed')")
    FailedJobs['RosCon2013_P16_score12prime'] = [r['ID'] for r in results]
    FailedExperiments1 = set([(r['ExperimentID'], r['UserDataSetExperimentID']) for r in results])

    results = ddGdb.execute("SELECT ID, ExperimentID, UserDataSetExperimentID FROM Prediction WHERE PredictionSet='RosettaCon2013_P16_talaris2013' AND Status IN ('failed')")
    FailedJobs['RosettaCon2013_P16_talaris2013'] = [r['ID'] for r in results]
    FailedExperiments2 = set([(r['ExperimentID'], r['UserDataSetExperimentID']) for r in results])

    results = ddGdb.execute("SELECT ID, ExperimentID, UserDataSetExperimentID FROM Prediction WHERE PredictionSet='RosCon2013_P16_talaris2013' AND Status IN ('failed')")
    FailedJobs['RosCon2013_P16_talaris2013'] = [r['ID'] for r in results]
    FailedExperiments3= set([(r['ExperimentID'], r['UserDataSetExperimentID']) for r in results])

    assert(FailedExperiments1 == FailedExperiments2)
    assert(FailedExperiments1.union(FailedExperiments3) == FailedExperiments1)

    print('FailedExperiments1 = %s' % str(list(FailedExperiments1)))
    print('FailedJobs = %s' % str(FailedJobs))

    check_dir = '/mnt/livewebserver/cluster/temp/'
    real_dir = '/var/cluster/temp/'
    import os
    for k, v in FailedJobs.iteritems():
        pth = os.path.join(check_dir, k)
        real_pth = os.path.join(real_dir, k)
        assert(os.path.exists(pth))
        for jobID in v:
            respath = os.path.join(pth, str(jobID))
            if os.path.exists(respath):
                respath = os.path.join(real_pth, str(jobID))
                print("sudo mv %s /var/cluster/temp/checkerrors/" % respath) # sudo -u klabqb3backrub
    sys.exit(0)


if False:
    for r in ddGdb.execute("SELECT ID, ExperimentID, PredictionSet, ddG FROM Prediction WHERE UserDataSetExperimentID=2698"):
        import klab.colortext as colortext
        colortext.message("%d, %s" % (r['ID'], r['PredictionSet']))
        import pickle
        ddG = pickle.loads(r['ddG'])
        for k, v in ddG['data'].iteritems():
            colortext.warning(k)
            print(v)

if False:
    print("***")
    for r in ddGdb.execute("SELECT ID, ExperimentID, PredictionSet, ddG FROM Prediction WHERE UserDataSetExperimentID=1588"):
        import klab.colortext as colortext
        colortext.message("%d, %s" % (r['ID'], r['PredictionSet']))
        import pickle
        ddG = pickle.loads(r['ddG'])
        for k, v in ddG['data'].iteritems():
            colortext.warning(k)
            print(v)

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

from klab import colortext
colortext.message("*** %s ***" % prediction_set_for_analysis)

score_cap = 10

if True:
	colortext.message("\n*** %s kellogg.total ***" % prediction_set_for_analysis)
	analyzer = analysis.Analyzer(prediction_set_for_analysis, ddG_score_type = 'kellogg.total', quiet_level = 1, score_cap = score_cap)
	analyzer.AddPublishedDDGsToAnalysisTables()
	reporter = analysis.Reporter(analyzer)
	reporter.CreateReport(description = analyzer.description, outfname = '%s_kellogg.pdf' % prediction_set_for_analysis, filetype = filetype)

if True:
#	for score_type in ['noah_6,0A', 'noah_7,0A', 'noah_8,0A', 'noah_9,0A']:
	for score_type in ['noah_8,0A']:
		#for score_method in ['total', 'positional', 'positional_twoscore']:
		for score_method in ['positional']:
			colortext.message("\n*** %s %s.%s ***" % (score_type, score_method, prediction_set_for_analysis))
			analyzer = analysis.Analyzer(prediction_set_for_analysis, ddG_score_type = '%s.%s' % (score_type, score_method), score_cap = score_cap)
			analyzer.AddPublishedDDGsToAnalysisTables()
			reporter = analysis.Reporter(analyzer)
			reporter.CreateReport(description = analyzer.description, outfname = '%s_%s_%s.pdf' % (prediction_set_for_analysis, score_type, score_method), filetype = filetype)

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
