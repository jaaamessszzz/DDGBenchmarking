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

    ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "RosettaCon2013_P16_talaris2013", "Protocol16 3.5.0 (talaris2013)", False,
            StoreOutput = True, Description = {}, InputFiles = {}, testonly = False,
            #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
    )

if False:

    ddG_connection.create_PredictionSet('Score12PrimeTest', halted = True, Priority = 10, BatchSize = 2)
    ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "Score12PrimeTest", "Protocol16 3.5.0 (score12prime)", False,
            StoreOutput = True, Description = {}, InputFiles = {}, testonly = False, shortrun = True
            #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
    )

if True:
    pass
    #ddG_connection.create_PredictionSet('RosettaCon2013_P16_score12prime', halted = True, Priority = 4, BatchSize = 40)
    #ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "RosettaCon2013_P16_score12prime", "Protocol16 3.5.0 (score12prime)", False,
    #        StoreOutput = True, Description = {}, InputFiles = {}, testonly = False, shortrun = False
    #        #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
    #)
    #ddG_connection.charge_PredictionSet_by_number_of_residues('RosettaCon2013_P16_score12prime')
    #Fix1TENEntries()
    #SetToActive()

if True:

    #ddG_connection.create_PredictionSet('RosCon2013_P16_score12prime', halted = True, Priority = 5, BatchSize = 2)
    #ddG_connection.create_PredictionSet('RosCon2013_P16_talaris2013', halted = True, Priority = 5, BatchSize = 2)

    #ddG_connection.create_PredictionSet('score12prime_test', halted = True, Priority = 4, BatchSize = 2)

    #ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "score12prime_test", "Protocol16 3.5.1 (score12prime)", False,
    #        StoreOutput = True, Description = {}, InputFiles = {}, testonly = False, shortrun = True
    #        #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
    #)



    #ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "RosCon2013_P16_score12prime", "Protocol16 3.5.1 (score12prime)", False,
    #        StoreOutput = True, Description = {}, InputFiles = {}, testonly = False, shortrun = False
    #        #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
    #)
    ddG_connection.charge_PredictionSet_by_number_of_residues('RosCon2013_P16_score12prime')

    #ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "RosCon2013_P16_talaris2013", "Protocol16 3.5.1 (talaris2013)", False,
    #        StoreOutput = True, Description = {}, InputFiles = {}, testonly = False, shortrun = False
    #        #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
    #)
    #ddG_connection.charge_PredictionSet_by_number_of_residues('RosCon2013_P16_talaris2013')

    #Fix num iterations to 50
    #Fix1TENEntries()
    #SetToActive()

