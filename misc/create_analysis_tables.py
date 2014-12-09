import sys
import os
#sys.path.insert(0, "..")
sys.path.insert(0, "../ddglib")
sys.path.insert(0, "../..")
sys.path.insert(0, "..")

from ddglib import analysis
from tools import colortext
from tools.fs.fsio import read_file, get_file_lines

filetype = "pdf"
prediction_set_for_analysis = 'Protocol_16_r57471'
prediction_set_for_analysis = 'RosCon2013_P16_score12prime'

colortext.message("*** %s ***" % prediction_set_for_analysis)

analyzer = analysis.Analyzer(prediction_set_for_analysis, ddG_score_type = 'kellogg.total', quiet_level = 1)
analyzer.AddPublishedDDGsToAnalysisTables()
for table_name, analysis_table in analyzer.analysis_tables.iteritems():
    if table_name not in []: # ProTherm, 'Potapov', 'Kellogg'
        colortext.message('Creating CSV for %s.' % table_name)
        csv_path = analyzer.CreateCSVFile(table_name, path = "/tmp")
        print(csv_path)
        #print(read_file(csv_path))
        c = 1
        record_numbers = []
        for l in get_file_lines(csv_path):
            print('%04d: %s' % (c, l))
            record_numbers.append(l.split(',')[5])
            c += 1
        print('%d lines.' % len(get_file_lines(csv_path)))
        print('%d lines, %d distinct record numbers.' % (len(record_numbers), len(set(record_numbers))))
        os.remove(csv_path)
        break