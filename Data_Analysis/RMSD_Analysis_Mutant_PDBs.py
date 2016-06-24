import pandas as pd
import sys
import numpy as np
import pprint
import seaborn as sns
import matplotlib.pyplot as plt
import pickle

def dataframe_construction(StructuralMetrics_pickle):
    mutant_df = pd.DataFrame(columns=('Prediction ID', 'WT PDBID', 'Mutant PDBID', 'RMSD Type', 'Point Mutant', 'RMSD'))
    PDB_set = set()
    PredID_set = set()
    mutation_list = []

    def add_to_df(PredictionID, Mutant_PDB, rmsd_type, point_mutant, rmsd_value, mutant_df):
        temp_df = pd.DataFrame(columns=('Prediction ID', 'WT PDBID', 'Mutant PDBID', 'RMSD Type', 'Point Mutant', 'RMSD'))
        temp_df.loc[Mutant_PDB.split()[0]] = pd.Series({'Prediction ID': PredictionID,
                                                        'WT PDBID': Mutant_PDB.split()[0],
                                                        'Mutant PDBID': Mutant_PDB,
                                                        'RMSD Type': rmsd_type[:-5],
                                                        'Point Mutant': point_mutant,
                                                        'RMSD': rmsd_value
                                                        }
                                                       )
        mutant_df = pd.concat([mutant_df, temp_df])
        return mutant_df

    with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/%s' %StructuralMetrics_pickle, 'rb') as input:
        RMSD_dict = pickle.load(input)

    Dataframe_pickle = '%s%s' % ('Dataframe-', StructuralMetrics_pickle[18:])

    try:
        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/%s' % Dataframe_pickle, 'rb') as input:
            mutant_df = pickle.load(input)
        for PredictionID in RMSD_dict:
            PredID_set.add(PredictionID)
            for nested_PredID in RMSD_dict[PredictionID]:
                for Mutant_PDB in RMSD_dict[PredictionID][nested_PredID]:
                    if RMSD_dict[PredictionID][nested_PredID][Mutant_PDB]['Global RMSD']['Mean'] > 15:
                        print '%s' % (Mutant_PDB)
                    else:
                        PDB_set.add(Mutant_PDB)

    except:
        for PredictionID in RMSD_dict:
            PredID_set.add(PredictionID)
            for nested_PredID in RMSD_dict[PredictionID]:
                for Mutant_PDB in RMSD_dict[PredictionID][nested_PredID]:
                    if RMSD_dict[PredictionID][nested_PredID][Mutant_PDB]['Global RMSD']['Mean'] > 15:
                        print '%s' %(Mutant_PDB)
                    else:
                        PDB_set.add(Mutant_PDB)
                        print Mutant_PDB
                        for rmsd_type in RMSD_dict[PredictionID][nested_PredID][Mutant_PDB]:
                            print rmsd_type
                            for rmsd in RMSD_dict[PredictionID][nested_PredID][Mutant_PDB][rmsd_type]:
                                if type(rmsd) == list: # rmsd[0] = Point Mutant Position, rmsd[1] = Dict
                                    mutation_list.append(rmsd[0])
                                    for rmsd_value in rmsd[1]['Raw']:
                                        mutant_df = add_to_df(PredictionID, Mutant_PDB, 'Point Mutant', rmsd[0], rmsd_value, mutant_df)
                                else:
                                    for rmsd_value in RMSD_dict[PredictionID][nested_PredID][Mutant_PDB][rmsd_type]['Raw']:
                                        mutant_df = add_to_df(PredictionID, Mutant_PDB, rmsd_type, None, rmsd_value, mutant_df)
        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/%s' % Dataframe_pickle, 'wb') as output:
            pickle.dump(mutant_df, output, 0)

    return RMSD_dict, mutant_df, PDB_set, mutation_list, PredID_set

def plot_stuff(mutant_df, PDB_set, mutation_list, PredID_set):

    # Import dataframes
    data_df = pd.read_csv('/kortemmelab/home/kyleb/reports/160608/analysis_sets/ZEMu/'
                         'ddg_analysis_type_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-score_method_Rescore-Talaris2014/data.csv')
    indexed_data_df = data_df.set_index('DatasetID')
    skempi_df = pd.read_csv('/kortemmelab/home/james.lucas/skempi_mutants.tsv', delimiter='\t')
    plotty_df = pd.DataFrame(columns=['Mutant PDBID', 'Predicted DDG', 'DDGs', 'RMSD Mean', 'RMSD StdDev'])

    df_agg = mutant_df.groupby('Mutant PDBID').agg({'RMSD': [np.mean, np.std]})
    df_agg.columns = ['RMSD Mean', 'RMSD StdDev']
    print df_agg

    # Modify skempi_df to include Predicted DDG values
    import re
    for index, row in skempi_df.iterrows():
        split_mutation = re.sub(':|-|>', ' ', row['Mutations']).split(',')
        temp_string_list = []
        for mutation in split_mutation:
            letters = mutation.split()
            temp_string_list.append('%s %s %s %s' %(letters[0], letters[1], ('   ' + letters[2])[-3:], letters[3]))
        joinme = '; '
        formatted_mutation = joinme.join(temp_string_list)
        for asdf, ghjk in data_df.iterrows():
            if formatted_mutation == ghjk['Mutations']:
                if row['Wildtype'] == ghjk['PDBFileID']:
                    temp_df = pd.Series([row['Mutant'], ghjk['Predicted'], row['DDGs'], df_agg.loc[(row['Wildtype'] + ' : ' + row['Mutant']), 'RMSD Mean'], df_agg.loc[row['Wildtype'] + ' : ' + row['Mutant'], 'RMSD StdDev']], index=['Mutant PDBID', 'Predicted DDG', 'DDGs', 'RMSD Mean', 'RMSD StdDev'])
                    plotty_df = plotty_df.append(temp_df, ignore_index = True)
                    print "It's a match!!!"
    print plotty_df

    for type, df_subset in mutant_df.groupby('RMSD Type'):
        if type != 'Point M':
            pass
            # fig, ax = plt.subplots()
            # ax = sns.boxplot(x='RMSD',y='Mutant PDBID',order=sorted(list(PDB_set)),data=df_subset)
            # plt.show()
            # fig.savefig("%s_boxplot_test.pdf" %type,)



    sys.exit(

    )
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = sns.jointplot(abs(plotty_df['DDGs'] - plotty_df['Predicted DDG']),
                       df_agg['RMSD Mean'],
                       xlim=(0, 8),
                       ylim=(0, 1.5),
                       size=10)
    plt.show()
    sys.exit()

    fig, ax = plt.subplots()
    g = sns.FacetGrid(mutant_df, col='WT PDBID', sharex=False, col_wrap=4, size=4, aspect=0.5)
    g.map(sns.boxplot, 'Mutant PDBID', 'RMSD')

    # ax = sns.boxplot(x=x,
    #                  y=y,
    #                  order=None,
    #                  data=df_subset)
    plt.show()
    # fig.savefig("%s_boxplot_test.pdf" %type,)


def main():
    StructuralMetrics_pickle = 'StructuralMetrics-ddg_analysis_type_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-score_method_Rescore-Talaris2014.pickle'
    RMSD_dict, mutant_df, PDB_set, mutation_list, PredID_set = dataframe_construction(StructuralMetrics_pickle)
    plot_stuff(mutant_df, PDB_set, mutation_list, PredID_set)


    # asdf_set = set()
    # for index, row in df.iterrows():
    #     asdf_set.add('%s : %s' %(row['Wildtype'], row['Mutant']))
    # print len(asdf_set)
    # print len(PDB_set)
    # print len(PredID_set)
    #
    # print asdf_set - PDB_set

main()