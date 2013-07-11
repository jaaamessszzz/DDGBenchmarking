import sys
import time
import profile

sys.path.insert(0, "..")
from tools import colortext
from tools.deprecated import rosettadb
from tools.debug.profile import ProfileTimer
from ddglib import dbapi, ddgdbapi
from ddglib import help as ddg_help
from ddglib.ddgfilters import *
from tools import pdb
import tools.deprecated.rosettahelper

ddG_connection = dbapi.ddG()

#Protocol16 3.5.0 (score12prime)
#Protocol16 3.5.0 (talaris2013)
#Protocol16 3.5.0 (baseline)
#Protocol16 3.5.0 (hbond_sp2_9g)
#Protocol16 3.5.0 (score12_hack_elec)

if False:
    ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "RosettaCon2013Protocol16_score12prime", "Protocol16 3.5.0 (score12prime)", False, StoreOutput = True, Description = {}, InputFiles = {}, testonly = False)
