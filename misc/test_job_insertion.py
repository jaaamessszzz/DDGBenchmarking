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
from ddglib.ppi_api import get_interface as get_ppi_interface
from ddglib.monomer_api import get_interface as get_protein_stability_interface


if __name__ == '__main__':

    # Set up the stability API
    stability_api = get_protein_stability_interface(read_file('ddgdb.pw'))

    # Set up the PPI API and get direct access to the database interface objects
    ppi_api = get_ppi_interface(read_file('ddgdb.pw'))
    ddg_db = ppi_api.DDG_db
    ddg_db_utf = ppi_api.DDG_db_utf

    # Informational_job tests
    ppi_api.get_prediction_set_details('RosCon2013_P16_score12prime')
    stability_api.get_prediction_ids('RosCon2013_P16_score12prime')

    # Print API help
    ppi_api.help()

    # Create the prediction set
    prediction_set_id = 'PPI test run'
    ppi_api.add_prediction_set(prediction_set_id, halted = True, priority = 7, batch_size = 41, allow_existing_prediction_set = True)
    ppi_api.alter_prediction_set_batch_size(prediction_set_id, 40)
    ppi_api.alter_prediction_set_priority(prediction_set_id, 5)

    pprint.pprint(stability_api.get_defined_user_datasets())
    pprint.pprint(ppi_api.get_defined_user_datasets())
    sys.exit(0)

    # Populate the prediction set with jobs
    ppi_api.add_prediction_run(prediction_set_id, 'AllBindingAffinityData', tagged_subset = 'ZEMu')

    for prediction_id in ppi_api.get_prediction_ids(prediction_set_id):
        print(prediction_id)
        # ... get job details, create command lines
        # ppi_api.add_job_command_lines(prediction_id, cmd_lines)

    sys.exit(0)
    a='''
    ppi_api.start_prediction_set(prediction_set_id)

    # Run the jobs
    while True:
        j = ppi_api.get_job(prediction_set)
        if j:
            ...create the jobs files on disk, group into an array job
        else:
            break
    #submit the job to the cluster, calling start_job(prediction_id) when the jobs have been submitted to SGE

    Inserting a job

    '''

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