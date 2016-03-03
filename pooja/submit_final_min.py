import os
import sys

#pdb_file=str(sys.argv[1])
#PDBID=pdb_file.split(".")

#json_mutation=sys.argv[2]

os.system("cp submit_minimize_last_wt.py wildtype/")
os.system("cp submit_minimize_last_mut.py mutation/")

os.chdir("wildtype/")
os.system("qsub submit_minimize_last_wt.py")
os.chdir("../mutation/")
os.system("qsub submit_minimize_last_mut.py")
os.chdir("../")
