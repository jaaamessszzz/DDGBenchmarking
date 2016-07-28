import pandas as pd
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import pickle
import re
import pprint
import math

def dataframe_construction(StructuralMetrics_pickle):
    mutant_df = pd.DataFrame(columns=('Prediction ID',
                                      'WT PDBID',
                                      'Mutant PDBID',
                                      'RMSD Type',
                                      'Point Mutant',
                                      'X1',
                                      'X2',
                                      'RMSD',
                                      'WT-Mutant Backbone RMSD',
                                      'Predicted DDG',
                                      'Experimental DDG',
                                      'Absolute Error DDG',
                                      'PDBID : Point Mutant',
                                      'Mutant Complex REU'
                                      )
                             )
    PDB_set = set()
    tossed_set = set()
    PredID_set = set()
    mutant_set = set()

    def add_to_df(PredictionID, Mutant_PDB, rmsd_type, point_mutant, rmsd_value, mutant_df, REU_score, backbone_rmsd, x1, x2):
        temp_df = pd.DataFrame(columns=('Prediction ID',
                                        'WT PDBID',
                                        'Mutant PDBID',
                                        'RMSD Type',
                                        'Point Mutant',
                                        'X1',
                                        'X2',
                                        'RMSD',
                                        'WT-Mutant Backbone RMSD',
                                        'PDBID : Point Mutant',
                                        'Mutant Complex REU'))

        temp_df.loc[Mutant_PDB.split()[0]] = pd.Series({'Prediction ID': PredictionID,
                                                        'WT PDBID': Mutant_PDB.split()[0],
                                                        'Mutant PDBID': Mutant_PDB,
                                                        'RMSD Type': rmsd_type.split()[2],
                                                        'Point Mutant': point_mutant,
                                                        'X1': x1,
                                                        'X2': x2,
                                                        'RMSD': rmsd_value,
                                                        'WT-Mutant Backbone RMSD': backbone_rmsd,
                                                        'PDBID : Point Mutant': '%s %s' %(Mutant_PDB, point_mutant),
                                                        'Mutant Complex REU': REU_score
                                                        }
                                                       )
        mutant_df = pd.concat([mutant_df, temp_df])

        return mutant_df

    with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/%s' %StructuralMetrics_pickle, 'rb') as input:
        RMSD_dict = pickle.load(input)

    Dataframe_pickle = '%s%s' % ('Dataframe-', StructuralMetrics_pickle[18:])

    try:
        print 'Opening %s\n' %Dataframe_pickle
        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/%s' % Dataframe_pickle, 'rb') as input:
            mutant_df = pickle.load(input)
        for PredictionID in RMSD_dict:
            PredID_set.add(PredictionID)
            if 'Failure' in RMSD_dict[PredictionID][PredictionID]:
                tossed_set.add(RMSD_dict[PredictionID][PredictionID])
            else:
                for Mutant_PDB in RMSD_dict[PredictionID][PredictionID]:
                    if RMSD_dict[PredictionID][PredictionID][Mutant_PDB]['Mutant - Global RMSD']['Mean'] > 10:
                        tossed_set.add(Mutant_PDB)
                    else:
                        PDB_set.add(Mutant_PDB)
                        for rmsd_type in RMSD_dict[PredictionID][PredictionID][Mutant_PDB]:
                            if 'Point_Mutant' in rmsd_type:
                                for point_mutant in RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type]:
                                    mutant_set.add('%s %s' % (Mutant_PDB, point_mutant))
    except:
        print '%s does not exist, generating it now...\n' %Dataframe_pickle
        for PredictionID in RMSD_dict:
            PredID_set.add(PredictionID)
            if 'Failure' in RMSD_dict[PredictionID][PredictionID]:
                tossed_set.add(RMSD_dict[PredictionID][PredictionID])
            else:
                for Mutant_PDB in RMSD_dict[PredictionID][PredictionID]:
                    bb_rmsd = RMSD_dict[PredictionID][PredictionID][Mutant_PDB]['WT-Mutant BackBone RMSDs']
                    if RMSD_dict[PredictionID][PredictionID][Mutant_PDB]['Mutant - Global RMSD']['Mean'] > 10:
                        tossed_set.add(Mutant_PDB)
                    else:
                        PDB_set.add(Mutant_PDB)
                        for rmsd_type in RMSD_dict[PredictionID][PredictionID][Mutant_PDB]:
                            if 'Point_Mutant' in rmsd_type:
                                for point_mutant in RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type]:
                                    mutant_set.add('%s %s' %(Mutant_PDB, point_mutant))
                                    for entry, rmsd_value in enumerate(RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type][point_mutant]['Raw']):
                                        if 'X1' in RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type][point_mutant]['X Angles']:
                                            if 'X2' in RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type][point_mutant]['X Angles']:
                                                mutant_df = add_to_df(PredictionID,
                                                                      Mutant_PDB,
                                                                      rmsd_type,
                                                                      point_mutant,
                                                                      rmsd_value,
                                                                      mutant_df,
                                                                      None,
                                                                      bb_rmsd['Point_Mutant'][point_mutant]['Mean'],
                                                                      RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type][point_mutant]['X Angles']['X1'][entry - 1],
                                                                      RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type][point_mutant]['X Angles']['X2'][entry - 1]
                                                                      )
                                            else:
                                                mutant_df = add_to_df(PredictionID,
                                                                      Mutant_PDB,
                                                                      rmsd_type,
                                                                      point_mutant,
                                                                      rmsd_value,
                                                                      mutant_df,
                                                                      None,
                                                                      bb_rmsd['Point_Mutant'][point_mutant]['Mean'],
                                                                      RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type][point_mutant]['X Angles']['X1'][entry - 1],
                                                                      None
                                                                      )
                                        else:
                                            mutant_df = add_to_df(PredictionID,
                                                                  Mutant_PDB,
                                                                  rmsd_type,
                                                                  point_mutant,
                                                                  rmsd_value,
                                                                  mutant_df,
                                                                  None,
                                                                  bb_rmsd['Point_Mutant'][point_mutant]['Mean'],
                                                                  None,
                                                                  None
                                                                  )


                            elif 'Global' in rmsd_type or 'Neighborhood' in rmsd_type:
                                for rmsd_value, REU_score in zip(RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type]['Raw'], RMSD_dict[PredictionID][PredictionID][Mutant_PDB]['Mutant Complex REUs']):
                                    mutant_df = add_to_df(PredictionID,
                                                          Mutant_PDB,
                                                          rmsd_type,
                                                          None,
                                                          rmsd_value,
                                                          mutant_df,
                                                          REU_score,
                                                          bb_rmsd[rmsd_type.split()[2]]['Mean'],
                                                          None,
                                                          None
                                                          )
                            else:
                                print rmsd_type
                                continue

        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/%s' % Dataframe_pickle, 'wb') as output:
            pickle.dump(mutant_df, output, 0)

    mutant_df = mutant_df.reset_index(drop = True)

    # Import dataframes
    # CHANGE FOR EACH RUN
    data_df = pd.read_csv('/kortemmelab/home/kyleb/reports/160608/analysis_sets/ZEMu/ddg_analysis_type_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-score_method_Rescore-Talaris2014/data.csv')
    # data_df = pd.read_csv('/kortemmelab/home/kyleb/reports/160627-bruball/analysis_sets/ZEMu/ddg_analysis_type_MatchPairs-prediction_set_id_zemu-psbrub_1.6-pv-nt50000-bruball-score_method_Rescore-Talaris2014/data.csv')
    skempi_df = pd.read_csv('/kortemmelab/home/james.lucas/skempi_mutants.tsv', delimiter='\t')

    # Modify mutant_df to include Predicted DDG values
    try:
        print 'Opening cached mutant_df\n'
        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/mutant_df_DDGs_cached-%s' %StructuralMetrics_pickle[18:], 'rb') as input:
            mutant_df = pickle.load(input)

    except:
        print 'Mutant_df cache does not exist, generating it now...\n'
        for skempi_index, skempi_row in skempi_df.iterrows():
            # Format mutation string
            split_mutation = re.sub(':|-|>', ' ', skempi_row['Mutations']).split(',')
            temp_string_list = []
            for mutation in split_mutation:
                letters = mutation.split()
                temp_string_list.append('%s %s %s %s' % (letters[0], letters[1], ('   ' + letters[2])[-3:], letters[3]))
            joinme = '; '
            formatted_mutation = joinme.join(temp_string_list)

            # Iterate through data_df to find correct DDG values
            for data_index, data_row in data_df.iterrows():
                if formatted_mutation == data_row['Mutations'] and skempi_row['Wildtype'] == data_row['PDBFileID'] and (
                        skempi_row['Wildtype'] + ' : ' + skempi_row['Mutant']) not in tossed_set:
                    # Add DDG values to mutant_df
                    print (skempi_row['Wildtype'] + ' : ' + skempi_row['Mutant'])
                    for mutant_index, mutant_row in mutant_df.iterrows():
                        if (skempi_row['Wildtype'] + ' : ' + skempi_row['Mutant']) == mutant_row['Mutant PDBID']:
                            mutant_df.loc[mutant_index, 'Predicted DDG'] = data_row[
                                                                               'Predicted'] / 1.2  # Scaling factor suggested by Shane
                            mutant_df.loc[mutant_index, 'Experimental DDG'] = data_row['Experimental_ZEMu']
                            mutant_df.loc[mutant_index, 'Absolute Error DDG'] = abs(
                                data_row['Experimental_ZEMu'] - data_row['Predicted'])
                    break

        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/mutant_df_DDGs_cached-%s' %StructuralMetrics_pickle[18:], 'wb') as output:
            pickle.dump(mutant_df, output, 0)

    # Generates bb_vs_sidechain_df dataframe
    RMSD_grouped = mutant_df.groupby('RMSD Type')
    neighborhood = RMSD_grouped.get_group('Neighborhood')
    point_mut = RMSD_grouped.get_group('Point_Mutant')

    try:
        print 'Opening cached bb_vs_sidechain_df\n'
        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/bb_vs_sidechain_df_cached-%s' % StructuralMetrics_pickle[18:], 'rb') as input:
            bb_vs_sidechain_df = pickle.load(input)

    except:
        print 'bb_vs_sidechain_df cache does not exist, generating it now...\n'
        bb_vs_sidechain_df = pd.DataFrame(columns=('Prediction ID',
                                                   'WT PDBID',
                                                   'Mutant PDBID',
                                                   'Point Mutant',
                                                   'Point Mutant RMSD',
                                                   'WT-Mutant Backbone RMSD',
                                                   'PDBID : Point Mutant',
                                                   'Experimental DDG',
                                                   'Predicted DDG',
                                                   'BB Group',
                                                   'DDG Group'
                                                   )
                                          )
        for name, group in point_mut.groupby('PDBID : Point Mutant'):
            for index, row in group.iterrows():
                if 'None' in row['PDBID : Point Mutant']:
                    pass
                else:
                    temp_series = pd.Series({'Prediction ID': row['Prediction ID'],
                                             'WT PDBID': row['WT PDBID'],
                                             'Mutant PDBID': row['Mutant PDBID'],
                                             'Point Mutant': row['Point Mutant'],
                                             'Point Mutant RMSD': row['RMSD'],
                                             'WT-Mutant Backbone RMSD':
                                                 neighborhood.groupby('Mutant PDBID').get_group(
                                                     row['Mutant PDBID'])['WT-Mutant Backbone RMSD'].iloc[0],
                                             'PDBID : Point Mutant': row['PDBID : Point Mutant'],
                                             'Experimental DDG': row['Experimental DDG'],
                                             'Predicted DDG': row['Predicted DDG']
                                             }
                                            )
                    bb_vs_sidechain_df = bb_vs_sidechain_df.append(temp_series, ignore_index=True)

        bb_vs_sidechain_df.sort_values('WT-Mutant Backbone RMSD', inplace=True)
        bb_vs_sidechain_df = bb_vs_sidechain_df.reset_index(drop=True)

        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/bb_vs_sidechain_df_cached-%s' % StructuralMetrics_pickle[18:], 'wb') as output:
            pickle.dump(bb_vs_sidechain_df, output, 0)

    # Assemble X1/X2 dataframes with neighborhood backbone RMSD column
    try:
        print 'Opening cached X angle df\n'
        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/X1_angle_df_cached-%s' % StructuralMetrics_pickle[18:], 'rb') as input:
            X1_angle_df = pickle.load(input)

    except:
        print 'X1_angle_df_cached does not exist, generating it now...\n'
        X1_angle_df = mutant_df[pd.notnull(mutant_df['X1'])]
        del X1_angle_df['WT-Mutant Backbone RMSD']  # These values are backbone RMSDs of single point mutant residues - we don't need that
        for index, row in X1_angle_df.iterrows():
            X1_angle_df.loc[index, 'BB RMSD'] = \
            neighborhood.groupby('Mutant PDBID').get_group(row['Mutant PDBID'])['WT-Mutant Backbone RMSD'].iloc[
                0]
            X1_angle_df.loc[index, 'Mutant X1'] = \
            RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']]['Mutant PDB X Angles'][
                row['Point Mutant']]['X1'][0]
            if 'X2' in RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']][
                'Mutant PDB X Angles'][row['Point Mutant']]:
                X1_angle_df.loc[index, 'Mutant X2'] = \
                RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']][
                    'Mutant PDB X Angles'][row['Point Mutant']]['X2'][0]

        X1_angle_df.sort_values('BB RMSD', inplace=True)
        X1_angle_df = X1_angle_df.reset_index(drop=True)

        with open('/kortemmelab/home/james.lucas/DDGBenchmarks_Test/Data_Analysis/RMSD_Outfiles/X1_angle_df_cached-%s' % StructuralMetrics_pickle[18:], 'wb') as output:
            pickle.dump(X1_angle_df, output, 0)


    # Output to csv for reference
    bb_vs_sidechain_df.to_csv('bb_vs_sidechain_df.csv')
    X1_angle_df.to_csv('X1_angle_df.csv')
    mutant_df.to_csv('Current_mutant_df.csv')

    return mutant_df, PDB_set, mutant_set, RMSD_dict, bb_vs_sidechain_df, X1_angle_df

def plot_stuff(mutant_df, PDB_set, mutant_set, RMSD_dict, bb_vs_sidechain_df, X1_angle_df):
    from matplotlib.backends.backend_pdf import PdfPages
    # output_pdf = PdfPages('RMSD_Analysis_Output_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000.pdf')
    # output_pdf = PdfPages('RMSD_Analysis_Output_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000_Current_Progress.pdf') # CHANGE FOR EACH RUN
    output_pdf = PdfPages('RMSD_Analysis_Output_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-X_Angle_Sorted.pdf') # CHANGE FOR EACH RUN

    print 'Unique WT PDBs: %s' %len(mutant_df['WT PDBID'].unique())
    print 'Unique Mutant PDBs: %s' %len(mutant_df['Mutant PDBID'].unique())
    print 'PDB_set length: %s' %len(PDB_set)

    ################################################
    # Backbone RMSD vs. Side Chain RMSD
    ################################################

    def BB_vs_Sidechain():
        # Make bins for BB RMSDs
        number_of_bins = 5
        bin_size = len(bb_vs_sidechain_df['WT-Mutant Backbone RMSD']) / number_of_bins + 1

        # Assign arbitrary bin identifiers for BB Group
        for index, row in bb_vs_sidechain_df.iterrows():
            bb_vs_sidechain_df.loc[index, 'BB Group'] = ((index + 1) // bin_size)
        # Find bin boundaries for BB group and add to dict
        bin_rename_dict = {}
        for name, group in bb_vs_sidechain_df.groupby('BB Group'):
            bin_rename_dict[name] = '%s -\n%s' % (group['WT-Mutant Backbone RMSD'].iloc[0], group['WT-Mutant Backbone RMSD'].iloc[len(group) - 1])
        # Rename bin identifiers to bin boundary values in BB group
        for index, row in bb_vs_sidechain_df.iterrows():
            bb_vs_sidechain_df.loc[index, 'BB Group'] = bin_rename_dict[bb_vs_sidechain_df.loc[index, 'BB Group']]

        # Assign bin identifiers for DDG Group
        for DDG_type in ['Experimental DDG', 'Predicted DDG']:
            for index, row in bb_vs_sidechain_df.iterrows():
                if row[DDG_type] > 2.5 or row[DDG_type] < -2.5:
                    bb_vs_sidechain_df.loc[index, DDG_type + ' Group'] = 'Extra Large DDG (DGG > 2.5 REU or DDG < -2.5 REU)'
                elif row[DDG_type] > 1 or row[DDG_type] < -1:
                    bb_vs_sidechain_df.loc[index, DDG_type + ' Group'] = 'Large DDG (2.5 REU > DGG > 1 REU or -2.5 < DDG < -1 REU)'
                elif row[DDG_type] > 0.5 or row[DDG_type] < -0.5:
                    bb_vs_sidechain_df.loc[index, DDG_type + ' Group'] = 'Medium DDG (1 REU > DGG > 0.5 REU or -1 < DDG < -0.5 REU)'
                else:
                    bb_vs_sidechain_df.loc[index, DDG_type + ' Group'] = 'Small DDG (0.5 REU > DDG > -0.5 REU)'

            sns.set_style('white', {'axes.grid': True, 'axes.edgecolor': '0'})
            sns.set_context('paper', font_scale=1.5, rc={'lines.linewidth': 1})

            fig, ax = plt.subplots(figsize=(20, 10))
            fig.suptitle('WT PDB - Mutant PDB Neighborhood Backbone RMSD vs. \nMutant PDB - RosettaOut Point Mutant Residues All-Atom RMSD', fontsize = 24, y=1.0)
            with sns.cubehelix_palette(number_of_bins, start=0.5, rot=-.75):
                sns.boxplot(x=bb_vs_sidechain_df['BB Group'],
                            y=bb_vs_sidechain_df['Point Mutant RMSD'],
                            ax=ax
                            )
            with sns.color_palette("husl", number_of_bins):
                sns.stripplot(x='BB Group',
                              y='Point Mutant RMSD',
                              hue= DDG_type + ' Group',
                              data=bb_vs_sidechain_df,
                              jitter=True,
                              ax=ax
                              )

            ax.set(xlabel='WT PDB - Mutant PDB Neighborhood Backbone RMSD', ylabel='Mutant PDB - RosettaOut Point Mutant Residues All-Atom RMSD')
            output_pdf.savefig(fig, pad_inches=1, bbox_inches='tight')

    ################################################
    # X angles
    ################################################
    def X_Angles_binned():
        # Make bins for BB RMSDs
        number_of_bins = 5
        X1_bin_size = len(X1_angle_df['BB RMSD']) / number_of_bins + 1

        # Assign arbitrary bin identifiers for BB Group
        for index, row in X1_angle_df.iterrows():
            X1_angle_df.loc[index, 'BB Group'] = ((index + 1) // X1_bin_size)

        # Find bin boundaries for BB group and add to dict
        bin_rename_dict = {}
        for name, group in X1_angle_df.groupby('BB Group'):
            bin_rename_dict[name] = '%s -\n%s' % (group['BB RMSD'].iloc[0], group['BB RMSD'].iloc[len(group) - 1])

        # Rename bin identifiers to bin boundary values in BB group
        for index, row in X1_angle_df.iterrows():
            X1_angle_df.loc[index, 'BB Group'] = bin_rename_dict[X1_angle_df.loc[index, 'BB Group']]

        # Categorize X1 angles as acceptable or not
        for index, row in X1_angle_df.iterrows():
            X1_angle = int(row['X1'])
            X1_upper_bound = int(row['Mutant X1'])+40
            X1_lower_bound = int(row['Mutant X1'])-40

            if X1_upper_bound > 180:
                if 180 >= X1_angle > X1_lower_bound or (-360 + X1_upper_bound) > X1_angle >= -180:
                    X1_angle_df.loc[index, 'X1 within 40'] = 1
                else:
                    X1_angle_df.loc[index, 'X1 within 40'] = 0
            elif X1_lower_bound < -180:
                if X1_upper_bound >= X1_angle > -180 or  180 > X1_angle >= (360 + X1_lower_bound):
                    X1_angle_df.loc[index, 'X1 within 40'] = 1
                else:
                    X1_angle_df.loc[index, 'X1 within 40'] = 0
            else:
                if X1_upper_bound > X1_angle > X1_lower_bound:
                    X1_angle_df.loc[index, 'X1 within 40'] = 1
                else:
                    X1_angle_df.loc[index, 'X1 within 40'] = 0

        # Categorize X1+X2 angles as acceptable or not
        for index, row in X1_angle_df.iterrows():
            if pd.isnull(row['X2']):
                pass
            else:
                X2_angle = int(row['X2'])
                X2_upper_bound = int(row['Mutant X2'])+40
                X2_lower_bound = int(row['Mutant X2'])-40

                if X2_upper_bound > 180:
                    if (180 >= X2_angle > X2_lower_bound or (-360 + X2_upper_bound) > X2_angle >= -180) and row['X1 within 40'] == 1:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 1
                    else:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 0
                elif X2_lower_bound < -180:
                    if (X2_upper_bound >= X2_angle > -180 or  180 > X2_angle >= (360 + X2_lower_bound)) and row['X1 within 40'] == 1:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 1
                    else:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 0
                else:
                    if X2_upper_bound > X2_angle > X2_lower_bound and row['X1 within 40'] == 1:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 1
                    else:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 0

        # Plotting!
        fig, ax = plt.subplots(figsize=(10, 10))
        sns.pointplot(x = 'BB Group',
                      y = 'X1 within 40',
                      data = X1_angle_df,
                      ci = None,
                      color = 'r',
                      ax=ax)
        sns.pointplot(x = 'BB Group',
                      y = 'X1+X2 within 40',
                      data = X1_angle_df,
                      ci = None,
                      color = 'b',
                      ax=ax)
        title = fig.suptitle('Point Mutant X-Angle Prediction',fontsize=24, y=.95)
        ax.set(xlabel='WT PDB - Mutant PDB Neighborhood Backbone RMSD',
               ylabel='Proportion of X angles Predicted within 40 Degrees')
        ax.set_ylim(0,1)

        # Proxy artists for legend...
        from matplotlib.lines import Line2D
        X1_handle = Line2D(range(1), range(1), color="red", marker='o', markerfacecolor="red")
        X1X2_handle = Line2D(range(1), range(1), color="blue", marker='o', markerfacecolor="blue")
        lgd = ax.legend([X1_handle, X1X2_handle], ['X1', 'X1 + X2'], loc='upper left')

        output_pdf.savefig(fig, pad_inches=1, bbox_extra_artists=[title, lgd], bbox_inches='tight')

    def X_Angles_sorted():
        # Categorize X1 angles as acceptable or not
        for index, row in X1_angle_df.iterrows():
            X1_angle = int(row['X1'])
            X1_upper_bound = int(row['Mutant X1'])+40
            X1_lower_bound = int(row['Mutant X1'])-40

            if X1_upper_bound > 180:
                if 180 >= X1_angle > X1_lower_bound or (-360 + X1_upper_bound) > X1_angle >= -180:
                    X1_angle_df.loc[index, 'X1 within 40'] = 1
                else:
                    X1_angle_df.loc[index, 'X1 within 40'] = 0
            elif X1_lower_bound < -180:
                if X1_upper_bound >= X1_angle > -180 or  180 > X1_angle >= (360 + X1_lower_bound):
                    X1_angle_df.loc[index, 'X1 within 40'] = 1
                else:
                    X1_angle_df.loc[index, 'X1 within 40'] = 0
            else:
                if X1_upper_bound > X1_angle > X1_lower_bound:
                    X1_angle_df.loc[index, 'X1 within 40'] = 1
                else:
                    X1_angle_df.loc[index, 'X1 within 40'] = 0

        # Categorize X1+X2 angles as acceptable or not
        for index, row in X1_angle_df.iterrows():
            if pd.isnull(row['X2']):
                pass
            else:
                X2_angle = int(row['X2'])
                X2_upper_bound = int(row['Mutant X2'])+40
                X2_lower_bound = int(row['Mutant X2'])-40

                if X2_upper_bound > 180:
                    if (180 >= X2_angle > X2_lower_bound or (-360 + X2_upper_bound) > X2_angle >= -180) and row['X1 within 40'] == 1:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 1
                    else:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 0
                elif X2_lower_bound < -180:
                    if (X2_upper_bound >= X2_angle > -180 or  180 > X2_angle >= (360 + X2_lower_bound)) and row['X1 within 40'] == 1:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 1
                    else:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 0
                else:
                    if X2_upper_bound > X2_angle > X2_lower_bound and row['X1 within 40'] == 1:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 1
                    else:
                        X1_angle_df.loc[index, 'X1+X2 within 40'] = 0

        X_set = set()
        for mutant_pdb in mutant_set:
            if mutant_pdb in X1_angle_df['PDBID : Point Mutant'].unique():
                X_set.add(mutant_pdb)

        X1_agg = X1_angle_df.groupby(['PDBID : Point Mutant']).aggregate(np.mean)
        X1_X2_agg = X1_agg[pd.notnull(X1_agg['X1+X2 within 40'])]

        X1_agg.to_csv('X1_agg.csv')

        # Grouping by WT PDB > Mutant PDB > Sorting by X1 scores
        for index, row in X1_agg.iterrows():
            X1_agg.loc[index, 'WT PDBID'] = index.split()[0]
            X1_agg.loc[index, 'Mutant PDBID'] = index.split()[2]
        WT_Mut_X1_sorted = []
        for WT_name, WT_subset in X1_agg.groupby(['WT PDBID']):
            for Mut_name, Mut_subset in WT_subset.groupby(['Mutant PDBID']):
                temp_order = sorted([(row['X1 within 40'], index) for index, row in Mut_subset.iterrows()])
                for asdf in temp_order:
                    WT_Mut_X1_sorted.append(asdf[1])

        print '\n\n\nWT_Mut_X1 SORTED!!!!!\n\n\n'
        pprint.pprint(WT_Mut_X1_sorted)
        print '\n\n\nWT_Mut_X1 SORTED!!!!!\n\n\n'

        # Sorting by X1 score
        jank_sort_list_X1 = sorted(list(set([(row['X1 within 40'], index) for index, row in X1_agg.iterrows()])))
        pprint.pprint(jank_sort_list_X1)
        jank_sort_list_X2 = sorted(list(set([(row['X1+X2 within 40'], index) for index, row in X1_X2_agg.iterrows()])))
        pprint.pprint(jank_sort_list_X2)

        pplot_X1_sort = [tup[1] for tup in jank_sort_list_X1]
        pprint.pprint(pplot_X1_sort)
        # Plotting!
        sns.set_style('white', {'axes.grid': True, 'axes.edgecolor': '0'})
        sns.set_context('poster', font_scale=1, rc={'lines.linewidth': 3})

        fig, ax = plt.subplots(figsize=(10, 10))
        sns.pointplot(x = 'PDBID : Point Mutant',
                      y = 'X1 within 40',
                      data = X1_angle_df,
                      palette='hls',
                      ci = None,
                      scatter_kws={"s": 200},
                      join=False,
                      order=pplot_X1_sort,
                      # order = sorted(list(X_set)),
                      # order=WT_Mut_X1_sorted,
                      hue='WT PDBID',
                      markers=['o']*len(WT_Mut_X1_sorted),
                      ax=ax
                      )
        sns.pointplot(x = 'PDBID : Point Mutant',
                      y = 'X1+X2 within 40',
                      data = X1_angle_df,
                      palette='hls',
                      ci = None,
                      scatter_kws={"s": 200},
                      join=False,
                      order=pplot_X1_sort,
                      # order=sorted(list(X_set)),
                      # order=WT_Mut_X1_sorted,
                      hue='WT PDBID',
                      markers=['*']*len(WT_Mut_X1_sorted),
                      ax=ax
                      )

        plt.xticks(rotation = 'vertical')

        title = fig.suptitle('Point Mutant X-Angle Prediction',fontsize=24, y=.95)
        ax.set(xlabel='WT PDB - Mutant PDB Neighborhood Backbone RMSD',
               ylabel='Proportion of X angles Predicted within 40 Degrees')
        ax.set_ylim(0,1)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend([handles[0], handles[len(handles)/2]],
                  ['X1', 'X1 + X2'],
                  loc='upper left',
                  fontsize=16,
                  )

        output_pdf.savefig(fig, pad_inches=0.25, bbox_extra_artists=[title], bbox_inches='tight')

    X_Angles_sorted()
    # X_Angles_binned()
    # BB_vs_Sidechain()
    # X_Angles()

# Generates all RMSD Type specific plots
    for type, df_subset in mutant_df.groupby('RMSD Type'):
        if type != 'Point_Mutant':

            # # DEVELOPMENT
            # continue

            if type == 'Global':
                description = 'Global C-alpha RMSDs'
            if type == 'Neighborhood':
                description = 'All-atom RMSDs for Residues within 8A of Mutation'

            ################################################
            # Pairplot
            ################################################
            def pairplot():
                sns.set_style('white', {'axes.grid': True, 'axes.edgecolor':'0'})
                sns.set_context('paper', font_scale=1.5, rc={'lines.linewidth': 1})
                # sns.despine()
                WT_pairplot = sns.pairplot(df_subset,
                                           vars = ['Experimental DDG', 'Predicted DDG', 'Absolute Error DDG', 'Mutant Complex REU', 'RMSD'],
                                           size =3,
                                           hue='WT PDBID',
                                           hue_order=sorted(list(mutant_df['WT PDBID'].unique())),
                                           kind='scatter',
                                           diag_kind='hist'
                                           )# .add_legend(bbox_to_anchor=(1.1, 0.5))
                lgd = WT_pairplot.fig.legend(handles=df_subset['WT PDBID'], labels=df_subset['WT PDBID'], bbox_to_anchor=(1.05, 0.5))
                output_pdf.attach_note('This pairplot compares the various numerical variables contained within the mutant_df dataframe. The following variables are compared in a pairwise fashion where hue is WT PDBID: Experimental DDG, Predicted DDG, Absolute Error DDG, Mutant Complex Rosetta Energy, and RMSD')
                title = WT_pairplot.fig.suptitle('PairPlot for %s' %description, fontsize =24, y=1.05)
                output_pdf.savefig(WT_pairplot.fig, pad_inches = 1, bbox_extra_artists = [title, lgd], bbox_inches='tight')

                Mut_pairplot = sns.pairplot(df_subset,
                                            vars = ['Experimental DDG', 'Predicted DDG', 'Absolute Error DDG', 'Mutant Complex REU', 'RMSD'],
                                            size =3,
                                            hue='Mutant PDBID',
                                            hue_order= sorted(list(mutant_df['Mutant PDBID'].unique())),
                                            kind='scatter',
                                            diag_kind='hist'
                                           )# .add_legend(bbox_to_anchor=(1.1, 0.5))
                lgd = Mut_pairplot.fig.legend(handles=df_subset['Mutant PDBID'], labels=df_subset['Mutant PDBID'], bbox_to_anchor=(1.05, 0.5))
                output_pdf.attach_note('This pairplot compares the various numerical variables contained within the mutant_df dataframe. The following variables are compared in a pairwise fashion where hue is Mutant PDBID: Experimental DDG, Predicted DDG, Absolute Error DDG, Mutant Complex Rosetta Energy, and RMSD')
                title = Mut_pairplot.fig.suptitle('PairPlot for %s' % description, fontsize=24, y=1.05)
                output_pdf.savefig(Mut_pairplot.fig, pad_inches = 1, bbox_extra_artists = [title, lgd], bbox_inches='tight')

            ################################################
            # Scatter plot - DDG Error vs. RMSD
            ################################################
            def DGG_Err_vs_RMSD():
                fig, ax = plt.subplots(figsize=(10, 10))
                ax = sns.regplot('Absolute Error DDG',
                                 'RMSD',
                                 data=df_subset,
                                 x_estimator=np.mean,
                                 fit_reg=False)

                ax.set_xlabel('DDG Absolute Error', fontsize=16)
                ax.set_ylabel('Average RMSD', fontsize=16)
                ax.set_xlim(left = 0)
                ax.set_ylim(bottom = 0)

                title = fig.suptitle(
                    'Average RMSD of a Rosetta-Generated Ensemble\nAgainst a Reference Mutant Crystal Structure',
                    # wrap=True,
                    fontsize=24,
                    y=1.00
                )
                ax.set_title(description,fontsize=18)
                output_pdf.attach_note('Error bars represent 95% confidence interval. Each point represents an ensemble for which we possess a reference mutant crystal structure')
                output_pdf.savefig(fig, pad_inches = 1, bbox_inches='tight')
                plt.close()

            ################################################
            # RMSD vs. REU Scatter plot  and RMSD Distributions
            ################################################
            def RMSD_vs_REU():
                for wt_pdb, wt_pdb_subset in df_subset.groupby('WT PDBID'):
                    sns.despine()
                    sns.set_style('white', {'axes.grid': True})
                    sns.set_context('notebook', font_scale=1, rc={'lines.linewidth': 1})
                    pairgrid = sns.PairGrid(vars = ['RMSD', 'Mutant Complex REU'],
                                            data=wt_pdb_subset,
                                            hue='Mutant PDBID',
                                            size = 6
                                            ).map_offdiag(plt.scatter).map_diag(sns.kdeplot).add_legend(bbox_to_anchor=(1.1, 0.5))
                    # pairgrid.set(xlim=(0, None))
                    title = pairgrid.fig.suptitle('%s vs.\nMutant Complex Rosetta Energy: %s' %(description, wt_pdb),
                                          fontsize=24,
                                          y=1.05)
                    lgd = pairgrid.fig.legend(handles=wt_pdb_subset['Mutant PDBID'], labels=wt_pdb_subset['Mutant PDBID'], bbox_to_anchor=(1.1, 0.5))
                    output_pdf.savefig(pairgrid.fig, pad_inches = 1, bbox_extra_artists = [lgd, title], bbox_inches='tight')
                    plt.close()

            ###############################################
            #Boxplot
            ###############################################
            def boxplots():
                fig, ax = plt.subplots(figsize=(20, 10))
                sns.boxplot(x='RMSD',
                            y='Mutant PDBID',
                            order=sorted(list(PDB_set)),
                            data=df_subset,
                            ax=ax)
                title = fig.suptitle(
                    'Distribution of %s for\nRosetta-Generated Ensembles vs. Reference Mutant Crystal Structure' % description,
                    fontsize=24,
                    y=1.00
                )
                sns.set_context('notebook', font_scale=1, rc={'lines.linewidth': 1})

                output_pdf.savefig(fig, pad_inches = 1, bbox_inches='tight')
                plt.close()

            # pairplot()
            # DGG_Err_vs_RMSD()
            # RMSD_vs_REU()
            # boxplots()

        else:
            ################################################
            # Binned Predicted DDG vs. Point Mutant RMSD
            ################################################
            def Predicted_DDG_vs_PointMut_RMSD():
                df_subset.sort_values('Predicted DDG', inplace=True)
                df_subset_index = df_subset.reset_index(drop=True)

                # Assign DDG bins in df_subset
                number_of_bins = 5
                DDG_bin_size = len(df_subset_index['Predicted DDG']) / number_of_bins + 1

                # Assign arbitrary bin identifiers for BB Group
                for index, row in df_subset_index.iterrows():
                    df_subset_index.loc[index, 'DDG Group'] = ((index + 1) // DDG_bin_size)

                # Find bin boundaries for BB group and add to dict
                bin_rename_dict = {}
                for name, group in df_subset_index.groupby('DDG Group'):
                    bin_rename_dict[name] = '%s --\n%s' % (group['Predicted DDG'].iloc[0], group['Predicted DDG'].iloc[len(group) - 1])

                # Rename bin identifiers to bin boundary values in BB group
                for index, row in df_subset_index.iterrows():
                    df_subset_index.loc[index, 'DDG Group'] = bin_rename_dict[df_subset_index.loc[index, 'DDG Group']]

                # Plotting!!!
                fig, ax = plt.subplots(figsize=(20, 10))
                with sns.color_palette("husl", number_of_bins):
                    sns.boxplot(x='DDG Group',
                                y='RMSD',
                                data=df_subset_index,
                                ax=ax)
                    sns.stripplot(x='DDG Group',
                                  y='RMSD',
                                  data=df_subset_index,
                                  jitter=True,
                                  color=sns.xkcd_rgb["charcoal grey"],
                                  ax=ax)

                # with sns.cubehelix_palette(number_of_bins, start=0.5, rot=-.75):
                #     sns.violinplot(x='DDG Group',
                #                    y='RMSD',
                #                    data=df_subset_index,
                #                    ax=ax)

                title = fig.suptitle('Rosetta-Generated Ensemble vs. Reference Mutant Crystal Structure\nPoint Mutant All-Atom RMSDs',
                                     fontsize=24,
                                     y=1.00)
                ax.set_xlabel('Predicted DDG', fontsize=16)
                ax.set_ylabel('Mutant PDB - RosettaOut Point Mutant All-Atom RMSD', fontsize=16)
                output_pdf.savefig(fig, pad_inches=1, bbox_inches='tight')

            ################################################
            # Scatter plot
            ################################################
            def scatter_PointMut():
                sns_fig = sns.lmplot('Absolute Error DDG',
                                     'RMSD',
                                     hue='Mutant PDBID',
                                     data=df_subset,
                                     x_estimator=np.mean,
                                     fit_reg=False,
                                     legend=False,
                                     # legend_out=True,
                                     size=10
                                     )
                sns_fig.fig.suptitle('All Mutant Resides All-Atom RMSD vs. Absolute Error in Predicted DDG for Each\n\
                Rosetta-Generated Ensemble/Reference Mutant Crystal Structure Pair')
                sns_fig = sns_fig.despine().set_axis_labels('Predicted DDG Absolute Error', 'Point Mutant All-atom RMSD')
                sns_fig.set(xlim=(0, None), ylim=(0, None))
                sns.set_style('white', {'axes.grid': True})
                sns.set_context('notebook', font_scale=1, rc={'lines.linewidth': 1})
                output_pdf.attach_note('Each point represents the average all-atom RMSD plus 95% confidence interval vs. the absolute error in predicted DDG for individual point mutant residues in the Rosetta-generated ensembles.')
                output_pdf.savefig(sns_fig.fig, pad_inches = 1, bbox_inches='tight')

            ################################################
            # Boxplot
            ################################################
            def boxplot_PointMut():
                for name, point_mutant_df in df_subset.groupby('WT PDBID'):
                    fig, ax = plt.subplots(figsize=(20, 10))

                    with sns.color_palette("husl"):
                        sns.boxplot(x='RMSD',
                                    y='PDBID : Point Mutant',
                                    order=sorted([asdf for asdf in mutant_set if name == asdf.split()[0]]),
                                    data=point_mutant_df,
                                    ax=ax)

                    sns.set_style('white', {'axes.grid': True})
                    sns.set_context('paper', font_scale=1.5, rc={'lines.linewidth': 2})

                    title = fig.suptitle(
                        'RMSDs for Individual Point Mutant Residues in %s\nRosetta-Generated Ensemble vs. Reference Mutant Crystal Structure' %name,
                        fontsize=24,
                        y=1.00)
                    ax.set_xlabel('RMSD', fontsize=16)
                    ax.set_ylabel('WT:Mut PDBID - Point Mutant', fontsize=16)

                    output_pdf.savefig(fig, pad_inches = 1, bbox_inches='tight')
                    plt.close()

            # Predicted_DDG_vs_PointMut_RMSD()
            # scatter_PointMut()
            # boxplot_PointMut()

    output_pdf.close()

def main():
    # CHANGE FOR EACH RUN
    StructuralMetrics_pickle = 'StructuralMetrics-ddg_analysis_type_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-score_method_Rescore-Talaris2014.pickle'
    # StructuralMetrics_pickle = 'StructuralMetrics-zemu-psbrub_1.6-pv-nt50000-bruball.pickle'

    mutant_df, PDB_set, mutant_set, RMSD_dict, bb_vs_sidechain_df, X1_angle_df = dataframe_construction(StructuralMetrics_pickle)
    plot_stuff(mutant_df, PDB_set, mutant_set, RMSD_dict, bb_vs_sidechain_df, X1_angle_df)

main()