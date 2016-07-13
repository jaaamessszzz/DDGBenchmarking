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
                            if type(RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type]) == list and len(
                                    RMSD_dict[PredictionID][PredictionID][Mutant_PDB][
                                        rmsd_type]) < 50:  # rmsd[0] = Point Mutant Position, rmsd[1] = Dict
                                for point_mutant in RMSD_dict[PredictionID][PredictionID][Mutant_PDB][rmsd_type]:
                                    mutant_set.add('%s %s' % (Mutant_PDB, point_mutant[0]))
                                    print '%s %s' % (Mutant_PDB, point_mutant[0])
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
                        print Mutant_PDB
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
    print 'Unique WT PDBs: %s' % len(mutant_df['WT PDBID'].unique())
    print 'Unique Mutant PDBs: %s' % len(mutant_df['Mutant PDBID'].unique())

    # Import dataframes
    # CHANGE FOR EACH RUN
    data_df = pd.read_csv('/kortemmelab/home/kyleb/reports/160608/analysis_sets/ZEMu/ddg_analysis_type_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-score_method_Rescore-Talaris2014/data.csv')
    # data_df = pd.read_csv('/kortemmelab/home/kyleb/reports/160627-bruball/analysis_sets/ZEMu/ddg_analysis_type_MatchPairs-prediction_set_id_zemu-psbrub_1.6-pv-nt50000-bruball-score_method_Rescore-Talaris2014/data.csv')
    skempi_df = pd.read_csv('/kortemmelab/home/james.lucas/skempi_mutants.tsv', delimiter='\t')

    # Modify mutant_df to include Predicted DDG values
    try:
        print 'Opening cached mutant_df\n'
        with open('mutant_df_DDGs_cached-%s' %StructuralMetrics_pickle[18:], 'rb') as input:
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

        with open('mutant_df_DDGs_cached-%s' %StructuralMetrics_pickle[18:], 'wb') as output:
            pickle.dump(mutant_df, output, 0)

    # # DEBUGGING
    import pprint
    pprint.pprint(mutant_df)
    pprint.pprint(RMSD_dict)
    mutant_df.to_csv('Current_mutant_df.csv')

    return RMSD_dict, mutant_df, PDB_set, mutant_set

def plot_stuff(mutant_df, PDB_set, mutant_set, RMSD_dict):
    from matplotlib.backends.backend_pdf import PdfPages
    # output_pdf = PdfPages('RMSD_Analysis_Output_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000.pdf')
    output_pdf = PdfPages('RMSD_Analysis_Output_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000_X-Angle_Tests.pdf') # CHANGE FOR EACH RUN

    print 'Unique WT PDBs: %s' %len(mutant_df['WT PDBID'].unique())
    print 'Unique Mutant PDBs: %s' %len(mutant_df['Mutant PDBID'].unique())
    print 'PDB_set length: %s' %len(PDB_set)

    ################################################
    # Backbone RMSD vs. Side Chain RMSD
    ################################################

    # Generates bb_vs_sidechain_df dataframe
    RMSD_grouped = mutant_df.groupby('RMSD Type')
    neighborhood = RMSD_grouped.get_group('Neighborhood')
    point_mut = RMSD_grouped.get_group('Point_Mutant')

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
                print row['PDBID : Point Mutant']
                print neighborhood.groupby('Mutant PDBID').get_group(row['Mutant PDBID'])['WT-Mutant Backbone RMSD'].iloc[0]
                print '\n'
                temp_series = pd.Series({'Prediction ID': row['Prediction ID'],
                                         'WT PDBID': row['WT PDBID'],
                                         'Mutant PDBID': row['Mutant PDBID'],
                                         'Point Mutant': row['Point Mutant'],
                                         'Point Mutant RMSD': row['RMSD'],
                                         'WT-Mutant Backbone RMSD': neighborhood.groupby('Mutant PDBID').get_group(row['Mutant PDBID'])['WT-Mutant Backbone RMSD'].iloc[0],
                                         'PDBID : Point Mutant': row['PDBID : Point Mutant'],
                                         'Experimental DDG': row['Experimental DDG'],
                                         'Predicted DDG': row['Predicted DDG']
                                         }
                                        )
                bb_vs_sidechain_df = bb_vs_sidechain_df.append(temp_series, ignore_index=True)

    bb_vs_sidechain_df.sort_values('WT-Mutant Backbone RMSD', inplace=True)
    bb_vs_sidechain_df = bb_vs_sidechain_df.reset_index(drop=True)

    # Output to csv for reference
    bb_vs_sidechain_df.to_csv('bb_vs_sidechain_df.csv')

    # Make bins for BB RMSDs
    number_of_bins = 5
    bin_size = len(bb_vs_sidechain_df['WT-Mutant Backbone RMSD']) / number_of_bins + 1

    print bin_size
    print len(bb_vs_sidechain_df['WT-Mutant Backbone RMSD'])

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
    for index, row in bb_vs_sidechain_df.iterrows():
        if row['Experimental DDG'] > 2.5 or row['Experimental DDG'] < -2.5:
            bb_vs_sidechain_df.loc[index, 'DDG Group'] = 'Extra Large DDG (DGG > 2.5 REU or DDG < -2.5 REU)'
        elif row['Experimental DDG'] > 1 or row['Experimental DDG'] < -1:
            bb_vs_sidechain_df.loc[index, 'DDG Group'] = 'Large DDG (2.5 REU > DGG > 1 REU or -2.5 < DDG < -1 REU)'
        elif row['Experimental DDG'] > 0.5 or row['Experimental DDG'] < -0.5:
            bb_vs_sidechain_df.loc[index, 'DDG Group'] = 'Medium DDG (1 REU > DGG > 0.5 REU or -1 < DDG < -0.5 REU)'
        else:
            bb_vs_sidechain_df.loc[index, 'DDG Group'] = 'Small DDG (0.5 REU > DDG > -0.5 REU)'

    # Plotting!
    sns.set_style('white', {'axes.grid': True, 'axes.edgecolor': '0'})
    sns.set_context('paper', font_scale=1.5, rc={'lines.linewidth': 1})

    fig, ax = plt.subplots(figsize=(20, 10))
    fig.suptitle('WT PDB - Mutant PDB Neighborhood Backbone RMSD vs. \nMutant PDB - RosettaOut Point Mutant Residues All-Atom RMSD', fontsize = 24)
    with sns.cubehelix_palette(number_of_bins, start=0.5, rot=-.75):
        sns.boxplot(x=bb_vs_sidechain_df['BB Group'],
                    y=bb_vs_sidechain_df['Point Mutant RMSD'],
                    ax=ax
                    )
    with sns.color_palette("husl", number_of_bins):
        sns.stripplot(x='BB Group',
                      y='Point Mutant RMSD',
                      hue='DDG Group',
                      data=bb_vs_sidechain_df,
                      jitter=True,
                      ax=ax
                      )

    ax.set(xlabel='WT PDB - Mutant PDB Neighborhood Backbone RMSD', ylabel='Mutant PDB - RosettaOut Point Mutant Residues All-Atom RMSD')

    output_pdf.savefig(fig, pad_inches=1, bbox_inches='tight')
    plt.show()

    # Scatterplot for Boxplot data
    box_scatter = sns.lmplot(x='WT-Mutant Backbone RMSD',
                             y='Point Mutant RMSD',
                             data = bb_vs_sidechain_df,
                             hue='DDG Group',
                             size = 10
                             )
    box_scatter.axes[0,0].set_ylim(0,)
    box_scatter.axes[0,0].set_xlim(0,)
    title = box_scatter.fig.suptitle('WT PDB - Mutant PDB Neighborhood Backbone RMSD vs. \nMutant PDB - RosettaOut Point Mutant Residues All-Atom RMSD', fontsize=24, y=1.05)
    output_pdf.savefig(box_scatter.fig, pad_inches=1, bbox_extra_artists=[title], bbox_inches='tight')
    plt.show()
    output_pdf.close()

    ################################################
    # X angles
    ################################################

    # Assemble X1/X2 dataframes with neighborhood backbone RMSD column
    X1_angle_df = mutant_df[pd.notnull(mutant_df['X1'])]
    del X1_angle_df['WT-Mutant Backbone RMSD'] # These values are backbone RMSDs of single point mutant residues - we don't need that
    X1_angle_df.to_csv('ashfjksdhfjksdhfkja.csv')
    for index, row in X1_angle_df.iterrows():
        X1_angle_df.loc[index, 'BB RMSD'] = neighborhood.groupby('Mutant PDBID').get_group(row['Mutant PDBID'])['WT-Mutant Backbone RMSD'].iloc[0]
        print RMSD_dict[row['Prediction ID']]
        print RMSD_dict[row['Prediction ID']][row['Prediction ID']]
        print RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']]
        print RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']]['Mutant PDB X Angles']
        print RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']]['Mutant PDB X Angles'][row['Point Mutant']]
        print RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']]['Mutant PDB X Angles'][row['Point Mutant']]['X1']
        X1_angle_df.loc[index, 'Mutant X1'] = RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']]['Mutant PDB X Angles'][row['Point Mutant']]['X1'][0]
        if 'X2' in RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']]['Mutant PDB X Angles'][row['Point Mutant']]:
            X1_angle_df.loc[index, 'Mutant X2'] = RMSD_dict[row['Prediction ID']][row['Prediction ID']][row['Mutant PDBID']]['Mutant PDB X Angles'][row['Point Mutant']]['X2'][0]

    X1_angle_df.to_csv('ashfjksdhfjksdhfkja.csv')
    sys.exit()

    # Make bins for BB RMSDs
    number_of_bins = 5
    X1_bin_size = len(X1_angle_df['BB RMSD']) / number_of_bins + 1

    print bin_size
    print len(bb_vs_sidechain_df['BB RMSD'])

    # Assign arbitrary bin identifiers for BB Group
    for index, row in X1_angle_df.iterrows():
        X1_angle_df.loc[index, 'BB Group'] = ((index + 1) // bin_size)
    # Find bin boundaries for BB group and add to dict
    bin_rename_dict = {}
    for name, group in X1_angle_df.groupby('BB Group'):
        bin_rename_dict[name] = '%s -\n%s' % (
        group['BB RMSD'].iloc[0], group['BB RMSD'].iloc[len(group) - 1])
    # Rename bin identifiers to bin boundary values in BB group
    for index, row in bb_vs_sidechain_df.iterrows():
        X1_angle_df.loc[index, 'BB Group'] = bin_rename_dict[X1_angle_df.loc[index, 'BB Group']]




    sys.exit()

    for type, df_subset in mutant_df.groupby('RMSD Type'):
        if type != 'Point_Mutant':

            if type == 'Global':
                description = 'Global C-alpha RMSDs'
            if type == 'Neighborhood':
                description = 'All-atom RMSDs for Residues within 8A of Mutation'

            continue

            ################################################
            # Pairplot
            ################################################
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
            for wt_pdb, wt_pdb_subset in df_subset.groupby('WT PDBID'):
                # gspec = gs.GridSpec(1, 2)
                # fig = plt.figure(figsize=(20, 10))
                # ax1 = fig.add_subplot(gspec[0,0])
                # ax2 = fig.add_subplot(gspec[0,1])
                #
                # sns.regplot('RMSD',
                #             'Mutant Complex REU',
                #             data=wt_pdb_subset,
                #             x_estimator=np.mean,
                #             scatter=True,
                #             fit_reg=False,
                #             ax=ax2
                #             )
                #
                # ax2.set_title(wt_pdb, fontsize=18)
                #
                # sns.distplot(wt_pdb_subset['RMSD'],
                #              hist=False,
                #              kde=True,
                #              norm_hist=True,
                #              ax=ax1)
                # ax1.set_ylabel('Percentage of cases (%)\nExcept not right now, working on it', fontsize=12)
                # ax1.set_xlim(left = 0)
                # ax1.set_ylim(bottom = 0)
                #
                # title = fig.suptitle(
                #     '%s vs.\nMutant Complex Rosetta Energy' %description,
                #     fontsize=24,
                #     y=1.00
                # )

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

        else:
            ################################################
            # X angles
            ################################################




            ################################################
            # Scatter plot
            ################################################
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
            for name, point_mutant_df in df_subset.groupby('WT PDBID'):
                fig, ax = plt.subplots(figsize=(20, 10))
                sns.set_style('white', {'axes.grid': True})
                sns.set_context('paper', font_scale=1.5, rc={'lines.linewidth': 2})

                sns.boxplot(x='RMSD',
                            y='PDBID : Point Mutant',
                            order=sorted([asdf for asdf in mutant_set if name == asdf.split()[0]]),
                            # hue='Mutant PDBID',
                            # hue_order=sorted(['%s : %s' %(asdf.split()[0], asdf.split()[2]) for asdf in mutant_set if name == asdf.split()[0]]),
                            data=point_mutant_df,
                            ax=ax)

                title = fig.suptitle(
                    'RMSDs for Individual Point Mutant Residues in %s\nRosetta-Generated Ensemble vs. Reference Mutant Crystal Structure' %name,
                    fontsize=24,
                    y=1.00)
                ax.set_xlabel('RMSD', fontsize=16)
                ax.set_ylabel('WT:Mut PDBID - Point Mutant', fontsize=16)

                output_pdf.savefig(fig, pad_inches = 1, bbox_inches='tight')
                plt.close()
    output_pdf.close()

def main():
    # CHANGE FOR EACH RUN
    StructuralMetrics_pickle = 'StructuralMetrics-ddg_analysis_type_CplxBoltzWT16.0-prediction_set_id_zemu-brub_1.6-nt10000-score_method_Rescore-Talaris2014.pickle'
    # StructuralMetrics_pickle = 'StructuralMetrics-zemu-psbrub_1.6-pv-nt50000-bruball.pickle'
    RMSD_dict, mutant_df, PDB_set, mutant_set = dataframe_construction(StructuralMetrics_pickle)
    plot_stuff(mutant_df, PDB_set, mutant_set, RMSD_dict)

main()