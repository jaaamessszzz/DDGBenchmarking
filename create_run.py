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
    #ddG_connection.charge_PredictionSet_by_number_of_residues('RosCon2013_P16_score12prime')

    #ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "RosCon2013_P16_talaris2013", "Protocol16 3.5.1 (talaris2013)", False,
    #        StoreOutput = True, Description = {}, InputFiles = {}, testonly = False, shortrun = False
    #        #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
    #)
    #ddG_connection.charge_PredictionSet_by_number_of_residues('RosCon2013_P16_talaris2013')

    #Fix num iterations to 50
    #Fix1TENEntries()
    #Postpone FailedExperiments
    #SetToActive()

    FailedExperiments = [(113528L, 4264L), (110786L, 1395L), (113529L, 4265L), (114163L, 4902L), (110325L, 897L), (114162L, 4901L), (110791L, 1400L), (112133L, 2776L), (110324L, 896L), (114150L, 4889L), (113515L, 4251L), (110789L, 1398L), (110321L, 893L), (113516L, 4252L), (114164L, 4903L), (114152L, 4891L), (114166L, 4905L), (110322L, 894L), (114151L, 4890L), (112131L, 2774L), (112454L, 3110L), (113532L, 4268L), (110317L, 889L), (112128L, 2771L), (113534L, 4270L), (112134L, 2777L), (109612L, 159L), (113536L, 4272L), (110790L, 1399L), (114160L, 4899L), (110580L, 1154L), (113525L, 4261L), (114161L, 4900L), (113689L, 4425L), (114155L, 4894L), (113523L, 4259L), (112630L, 3299L), (110319L, 891L), (113688L, 4424L), (112127L, 2770L), (113522L, 4258L), (113533L, 4269L), (110787L, 1396L), (113530L, 4266L), (113531L, 4267L), (113521L, 4257L), (114158L, 4897L), (110320L, 892L), (114157L, 4896L), (110318L, 890L), (112629L, 3298L), (113526L, 4262L), (110788L, 1397L), (113537L, 4273L), (110579L, 1153L), (113527L, 4263L), (112129L, 2772L), (110792L, 1401L), (110323L, 895L), (112135L, 2778L), (114159L, 4898L), (113518L, 4254L), (113538L, 4274L), (114154L, 4893L), (113519L, 4255L), (114165L, 4904L), (113513L, 4249L), (112132L, 2775L), (113520L, 4256L), (112455L, 3111L), (114156L, 4895L), (113514L, 4250L), (110578L, 1151L), (113535L, 4271L), (113524L, 4260L), (112130L, 2773L), (109613L, 160L), (113517L, 4253L)]
    assert(len(FailedExperiments) == 78)

    if False:
        ddG_connection.create_PredictionSet('RosCon2013_P16_talaris2013sc_test', halted = True, Priority = 5, BatchSize = 10)
        ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "RosCon2013_P16_talaris2013sc_test", "Protocol16 3.5.1 (talaris2013sc)", False,
                StoreOutput = True, Description = {}, InputFiles = {}, testonly = False, shortrun = True
                #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
        )
        ddG_connection.ddGDB.execute("UPDATE PredictionSet SET Status='active' WHERE ID='RosCon2013_P16_talaris2013sc_test' ")

    if False:
        ddG_connection.create_PredictionSet('RosCon2013_P16_talaris2013sc', halted = True, Priority = 5, BatchSize = 10)
        ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "RosCon2013_P16_talaris2013sc", "Protocol16 3.5.1 (talaris2013sc)", False,
                StoreOutput = True, Description = {}, InputFiles = {}, testonly = False, shortrun = False
                #only_single_mutations = False, # we do not need to run multiple mutations for RosettaCon 2013
        )
        ddG_connection.charge_PredictionSet_by_number_of_residues('RosCon2013_P16_talaris2013sc')
        # postpone experiments that will fail
        for fe in FailedExperiments:
            assert(len(ddG_connection.ddGDB.execute("SELECT ID FROM Prediction WHERE PredictionSet='RosCon2013_P16_talaris2013sc' AND ExperimentID=%s AND UserDataSetExperimentID=%s",
                                         parameters=(fe[0], fe[1]))) == 1)
            ddG_connection.ddGDB.execute("UPDATE Prediction SET Status='postponed' WHERE PredictionSet='RosCon2013_P16_talaris2013sc' AND ExperimentID=%s AND UserDataSetExperimentID=%s",
                                         parameters=(fe[0], fe[1]))

            assert(len(ddG_connection.ddGDB.execute("SELECT ID FROM Prediction WHERE PredictionSet='RosCon2013_P16_talaris2013' AND ExperimentID=%s AND UserDataSetExperimentID=%s",
                                         parameters=(fe[0], fe[1]))) == 1)
            ddG_connection.ddGDB.execute("UPDATE Prediction SET Status='postponed' WHERE PredictionSet='RosCon2013_P16_talaris2013' AND ExperimentID=%s AND UserDataSetExperimentID=%s",
                                         parameters=(fe[0], fe[1]))

        #Fix1TENEntries()
        #SetToActive()


    if False:
        results = ddG_connection.ddGDB.execute("SELECT ID, ExperimentID, ddG from Prediction WHERE PredictionSet='RosCon2013_P16_score12prime' AND Status='done'")
        results = ddG_connection.ddGDB.execute("SELECT ID, ExperimentID, ddG from Prediction WHERE ID=36237")
        results = ddG_connection.ddGDB.execute("SELECT ID, ExperimentID, ddG, PredictionSet from Prediction WHERE ExperimentID=112055 AND Status='done'")
        for r in results:
            ddGo = pickle.loads(r['ddG'])
            ddG = ddGo['data']['kellogg']['total']['ddG']
            print("")
            colortext.message(r['ID'])
            if ddG > 100 or True:
                ExperimentID = r['ExperimentID']
                PredictionSet = r['PredictionSet']
                PDBFileID = ddG_connection.ddGDB.execute("SELECT PDBFileID FROM Experiment WHERE ID=%s", parameters=(ExperimentID,))[0]['PDBFileID']
                print('Outlier: Prediction #%d, PredictionSet %s, ddG=%0.2f, ExperimentID=%d, PDBFileID=%s.' % (r['ID'], PredictionSet, ddG, ExperimentID, PDBFileID))
                mutations = ddG_connection.ddGDB.execute("SELECT * FROM ExperimentMutation WHERE ExperimentID=%s", parameters=(ExperimentID,))
                for m in mutations:
                    print('Chain %s: wildtype aa %s, residue ID %s, mutant aa %s ' % (m['Chain'], m['WildTypeAA'], m['ResidueID'], m['MutantAA']))


