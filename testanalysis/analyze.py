import sys
sys.path.insert(0, "..")
sys.path.insert(0, "../common")
sys.path.insert(0, "../ddglib")
from common import colortext, rosettadb
import analysis

analyzer = analysis.Analyzer("AllExperimentsProtocol16")
analyzer.AddPublishedDDGsToAnalysisTables()
reporter = analysis.Reporter(analyzer)
reporter.CreateReport(description = analyzer.description)


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
