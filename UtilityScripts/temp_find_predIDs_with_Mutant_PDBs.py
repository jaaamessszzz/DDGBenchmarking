import sys
import os
import pandas as pd

def main():
    # List of dictionaries for Resfile mutation info
    sys.path.insert(0,'/kortemmelab/home/james.lucas')  # this should point to the directory above your ddg repo checkout
    from kddg.api.ppi import get_interface_with_config_file as get_ppi_interface

    # Create an interface to the database
    ppi_api = get_ppi_interface()

    # Grabs Prediction IDs from .csv file
    csv_info = pd.read_csv('/kortemmelab/home/james.lucas/zemu-psbrub_1.6-pv-1000-lbfgs_IDs.csv')
    PredID_list = []
    for index, row in csv_info.iterrows():
        PredID_list.append(int(row[0].split()[0]))

    print PredID_list
    sys.exit()

    # Get details back for one prediction
    for predID in PredID_list
        PredID_Details = ppi_api.get_job_details(predID, include_files=False)

    # Something
    df = pd.read_csv('/kortemmelab/home/james.lucas/skempi_mutants.tsv', delimiter='\t')
    for index, row in df.iterrows():
        if row['PPMutagenesisID'] == PredID_Details['PDBMutations'][0]['PPMutagenesisID']:
            Mutant_PDB_ID = row['Mutant']
            break