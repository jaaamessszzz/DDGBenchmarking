import sys
import traceback
import datetime
import string
import os
import pickle
import numpy
import json
import pprint

if __name__ == "__main__":
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, "../updatedb")

import tools.colortext as colortext
from tools.fs.fsio import read_file
from ddglib.ppi_api import BindingAffinityDDGInterface
from ddglib.ppi_api import get_interface as get_ppi_interface
from ddglib.monomer_api import MonomericStabilityDDGInterface
from ddglib.monomer_api import get_interface as get_protein_stability_interface

if __name__ == '__main__':
    from ddglib.ppi_api import get_interface as get_ppi_interface
    ppi_api = get_ppi_interface(read_file('ddgdb.pw'))
    ppi_api.help()
    sys.exit(0)
    stability_api = get_protein_stability_interface(read_file('ddgdb.pw'))
    pprint.pprint(ppi_api.get_prediction_set_details('RosCon2013_P16_score12prime'))
    print(stability_api.get_prediction_ids('RosCon2013_P16_score12prime'))

    #from ddglib.monomer_api import get_interface as get_protein_stability_interface
    #stability_api = get_protein_stability_interface(read_file('ddgdb.pw'))
    #stability_api.help()

    #ddg_db = ddg_api.DDG_db
    #ddg_db_utf = ddg_api.DDG_db_utf

    a='''
    1.
    ppi_api.create_prediction_set("Shane's test prediction set")

    2.
    ppi_api.create_predictions_from_userdataset("Shane's test prediction set", 'AllBindingAffinity', tagged_subset = 'ZEMu')

    # Create a list of PredictionPPI records with:
    - PDB file with stripped chains
       - PPMutagenesis.ID
       - UserPPDataSetExperimentID (specifies PPMutagenesisID and PDB complex definition (PDB ID, PPComplexID, SetNumber))
       - ProtocolID none at present
       - Cost (num residues in stripped PDB)
       - KeptHETATMLines?
       - ResidueMapping (JSON from Rosetta numbering to PDB numbering)
       - InputFiles - mutfile/resfile?
       - Description
       - ScoreVersion
       - ddG (NULL from 23505 to 76632) - 1860 records
       - Scores (NULL from 23505 to 76632)
       - StructureScores (only non-NULL on records 55808, 55809 )

       Each prediction has a set of PredictionStructureScores
         - per prediction, score method (Global p16, Local 8A Noah, ...), score type (DDG, Mutant, Wildtype), run number e.g. 1-50, we store:
           - score components
           - DDG

    3.
    Kyle's runner populates these records with command lines?
    Kyle's runner runs the jobs and saves the results back into the database


    '''