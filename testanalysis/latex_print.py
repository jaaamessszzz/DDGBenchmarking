import sys
sys.path.insert(0, "../..")
import tools.latex as latex

code = '''
// blah blah blah
'''

codeprinter = latex.LaTeXCodePrinter('rosetta\_backend', code)
codeprinter.compile_pdf("codetest.pdf")