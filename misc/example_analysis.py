import sys
import os
import pprint
#sys.path.insert(0, "..")
sys.path.insert(0, "../ddglib")
sys.path.insert(0, "../..")
sys.path.insert(0, "..")

from ddglib import analysis
import pprint
from tools import colortext
from tools.fs.fsio import read_file, get_file_lines, write_file
import json

filetype = "pdf"

from ddglib import db_api, ddgdbapi, monomer_api
from tools.bio import pdb



def dump_data(prediction_set, outfile):
    ddG_connection = db_api.ddG()
    ddGdb = ddgdbapi.ddGDatabase()

    userdata_set = 'AllValidPGPK'

    cached_pdb_details = json.loads(read_file('cached_pdb_details.json'))
    analysis_breakdown = ddG_connection.get_predictionset_data(prediction_set, userdata_set, cached_pdb_details = cached_pdb_details)

    test_data = dict(
        amino_acids = analysis_breakdown.amino_acids,
        pdb_details = analysis_breakdown.pdb_details,
        predictions = analysis_breakdown.predictions,
        #single_mutation_GP_predictions = analysis_breakdown.single_mutation_GP_predictions,
        #single_mutation_no_GP_predictions = analysis_breakdown.single_mutation_no_GP_predictions,
        #multiple_mutation_predictions = analysis_breakdown.multiple_mutation_predictions,
        analysis_datasets = analysis_breakdown.analysis_datasets,
    )
    write_file(outfile, json.dumps(test_data))


def load_data(infile):
    d = json.loads(read_file(infile))
    return monomer_api.AnalysisBreakdown(d['amino_acids'], d['pdb_details'], d['predictions'], d['analysis_datasets'])

def get_data_for_small_large_diagram_for_website():
    d = json.loads(read_file('r57471_analysis_input.json'))

    ddG_connection = db_api.ddG()
    amino_acids = ddG_connection.get_amino_acids_for_analysis()

    amino_acid_volumes = {}
    for aa, details in amino_acids.iteritems():
        amino_acid_volumes[aa] = details['van_der_Waals_volume']

    assert(len(amino_acid_volumes) == 20)

    ddGdb = ddgdbapi.ddGDatabase()
    datasets = ['CuratedProTherm_2014/12/04', 'Guerois_10.1016/S0022-2836(02)00442-4_2002/07/05', 'Kellogg_10.1002/prot.22921_2010/12/03', 'Potapov_10.1093/protein/gzp030_2009/09/01']
    multiple_mutations = dict.fromkeys(datasets, 0)
    net_counts = dict.fromkeys(datasets, 0)
    SL_counts = {}
    for dataset in datasets:
        SL_counts[dataset] = {'SL': 0, 'LS': 0, 'XX': 0}

    for dataset in datasets:
        records = ddGdb.execute_select('SELECT * FROM DataSetDDG WHERE DataSetID=%s', parameters=(dataset,))
        print('%d records in %s' % (len(records), dataset))
        for r in records:
            experiment_ids = set([s['ExperimentID'] for s in ddGdb.execute_select('''
                SELECT ExperimentID FROM
                DataSetDDGSource
                INNER JOIN ExperimentAssayDDG ON DataSetDDGSource.ExperimentAssayID=ExperimentAssayDDG.ExperimentAssayID AND DataSetDDGSource.Type=ExperimentAssayDDG.Type
                INNER JOIN ExperimentAssay ON ExperimentAssayDDG.ExperimentAssayID=ExperimentAssay.ID
                WHERE DataSetDDGID=%s''', parameters = (r['ID'],))])
            if not len(experiment_ids) == 1:
                colortext.warning('Duplicate record in %s: Dataset record #%d, ExperimentIDs=%s.' % (dataset, r['ID'], ', '.join(map(str, sorted(experiment_ids)))))
                continue
            experiment_id = experiment_ids.pop()
            mutations = ddGdb.execute_select('''SELECT * FROM ExperimentMutation WHERE ExperimentID=%s''', parameters=(experiment_id,))
            if len(mutations) > 1:
                mutation_classes = set()
                error = False
                for mutation in mutations:
                    wt, mut =  mutation['WildTypeAA'], mutation['MutantAA']
                    if r['MutationIsReversed']:
                        # Note: For reverse mutations, we need to switch the order since we only store the forward mutation
                        wt, mut =  mutation['MutantAA'], mutation['WildTypeAA']

                    if wt == mut:
                        colortext.warning('Error in %s: Record mutating %s to %s in Experiment #%d.' % (dataset, wt, mut, experiment_id))
                        error = True
                    elif amino_acid_volumes[wt] < amino_acid_volumes[mut]:
                        mutation_classes.add('SL')
                    elif amino_acid_volumes[wt] > amino_acid_volumes[mut]:
                        mutation_classes.add('LS')
                    else:
                        assert(amino_acid_volumes[wt] == amino_acid_volumes[mut])
                        mutation_classes.add('XX')
                if not(error) and len(mutation_classes) == 1:
                    colortext.printf('Multiple mutation case allowed since both mutations have the same class.', 'cyan')
                    SL_counts[dataset][mutation_classes.pop()] += 1
                else:
                    multiple_mutations[dataset] += 1
                    continue # skip multiple mutations
            else:
                assert(len(mutations) == 1)
                mutation = mutations[0]

                wt, mut =  mutation['WildTypeAA'], mutation['MutantAA']
                if r['MutationIsReversed']:
                    # Note: For reverse mutations, we need to switch the order since we only store the forward mutation
                    wt, mut =  mutation['MutantAA'], mutation['WildTypeAA']

                if wt == mut:
                    colortext.warning('Error in %s: Record mutating %s to %s in Experiment #%d.' % (dataset, wt, mut, experiment_id))
                    continue
                elif amino_acid_volumes[wt] < amino_acid_volumes[mut]:
                    SL_counts[dataset]['SL'] += 1
                elif amino_acid_volumes[wt] > amino_acid_volumes[mut]:
                    SL_counts[dataset]['LS'] += 1
                else:
                    assert(amino_acid_volumes[wt] == amino_acid_volumes[mut])
                    SL_counts[dataset]['XX'] += 1
            net_counts[dataset] += 1
    #GASCPDTNVEQHMLIKFYRW
    colortext.message('\nRecords with multiple mutations that were skipped.')
    pprint.pprint(multiple_mutations)

    colortext.message('\nNet SL, LS, and XX counts for the datasets.')
    pprint.pprint(SL_counts)
    for dataset, details in SL_counts.iteritems():
        for type, type_total in details.iteritems():
            details[type] = 100 * (float(type_total)/float(net_counts[dataset]))
    colortext.message('\nNet SL, LS, and XX percentages for the datasets.')
    pprint.pprint(SL_counts)



def create_Venn_diagrams_for_website():
    import matplotlib
    from matplotlib import pyplot as plt
    import numpy as np
    from matplotlib_venn import venn3, venn3_circles

    figure, axes = plt.subplots(1, 2)
    #fig = plt.figure(figsize=(6,6))

    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 22}

    matplotlib.rc('font', **font)


    v1 = venn3(subsets=(1061, 8, 14, 141, 102, 12, 702), set_labels = ('Potapov', 'Guerois', 'Kellogg'), set_colors=('red', 'yellow', 'blue'), ax=axes[0])

    v2 = venn3(subsets=(1, 1, 88, 76, 1, 1, 136), set_labels = ('Potapov', 'Guerois', 'Kellogg'), set_colors=('red', 'yellow', 'blue'), ax=axes[1])
    v2.get_label_by_id('001').set_alpha(1.0)
    v2.get_label_by_id('001').set_color('blue')
    v2.get_label_by_id('010').set_alpha(1.0)
    v2.get_label_by_id('010').set_color('yellow')
    v2.get_label_by_id('011').set_alpha(1.0)
    v2.get_label_by_id('100').set_alpha(1.0)
    v2.get_label_by_id('100').set_color('red')
    v2.get_label_by_id('101').set_alpha(1.0)
    v2.get_label_by_id('110').set_alpha(1.0)
    v2.get_label_by_id('111').set_alpha(1.0)
    v2.get_label_by_id('100').set_text('')
    v2.get_label_by_id('011').set_text('')
    v2.get_label_by_id('101').set_text('')

    import matplotlib.patches as mpatches
    potapov_patch = mpatches.Patch(color=(1, 153.0/255.0, 153.0/255.0), label='Potapov')
    guerois_patch = mpatches.Patch(color=(1, 1, 153.0/255.0), label='Guerois')
    kellogg_patch = mpatches.Patch(color=(153.0/255.0, 153.0/255.0, 1), label='Kellogg')
    plt.legend(handles=[potapov_patch, kellogg_patch, guerois_patch], loc = (-0.3,-0.1), bbox_transform = plt.gcf().transFigure)

    #plt.legend( [1,2,3], ['A','B','C'], loc = 'lower center', bbox_to_anchor = (0,-0.1,1,1),
    #        bbox_transform = plt.gcf().transFigure )

    #plt.title("Intersection between the prior datasets and ProTherm*")
    plt.annotate('ProTherm*', xy=v2.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-100,-30), fontsize = 26,
                 ha='center', textcoords='offset points')
    plt.show()

    figure = plt.gcf()
    figure.set_size_inches(18.5,10.5)

    figure.savefig('test.png', dpi=100)

    return
    fig = plt.figure(figsize=(6,6))
    #v = venn3(subsets=(1061, 8, 14, 141, 102, 12, 702), set_labels = ('Potapov', 'Guerois', 'Kellogg'))
    #v = venn3(subsets=(0, 1, 88, 76, 0, 0, 136), set_labels = ('Potapov', 'Guerois', 'Kellogg'))
    #v = venn3(subsets=(1, 1, 88, 76, 1, 1, 136), set_labels = ('Potapov', 'Guerois', 'Kellogg'), set_colors=('blue', 'lightblue'))
    v = venn3(subsets=(1, 1, 88, 76, 1, 1, 136), set_labels = ('Potapov', 'Guerois', 'Kellogg'), set_colors=('red', 'yellow', 'blue'))
    #v.get_patch_by_id('100').set_alpha(1.0)
    #v.get_patch_by_id('100').set_color('white')
    v.get_label_by_id('100').set_text('')
    v.get_label_by_id('011').set_text('')
    v.get_label_by_id('101').set_text('')

    v.get_patch_by_id('110').set_color('red')

    #c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='solid', linewidth=1.0)
    #c[0].set_lw(1.0)
    #c[0].set_ls('solid')
    plt.title("Intersection between the prior datasets and ProTherm*")
    #plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
    #             ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    #             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    plt.show()

    fig.savefig('test.png')

if __name__ == '__main__':

    #get_data_for_small_large_diagram_for_website()
    #sys.exit(0)

    if False:
        import json
        print('Dumping data...')
        dump_data('RosCon2013_P16_score12prime', 'score12_analysis_input.json')
        dump_data('RosCon2013_P16_talaris2013', 'r55534_talaris2013_analysis_input.json')
        dump_data('RosCon2013_P16_talaris2013sc', 'r55534_talaris2013sc_analysis_input.json')
        dump_data('Protocol_16_r57471', 'r57471_analysis_input.json')
        print('Done')
        sys.exit(0)

    if False:
        # code used to create the benchmark capture
        ddG_connection = db_api.ddG()
        ddGdb = ddgdbapi.ddGDatabase()

        print(ddG_connection.determine_best_pair(72856, ScoreMethodID = 1))

        scores = {}
        wts = 0
        muts = 0
        records = ddGdb.execute_select('SELECT StructureID, total, ScoreType FROM `PredictionStructureScore` WHERE `PredictionID`=72856')
        for r in records:
            if r['ScoreType'] == 'WildType':
                scores[r['StructureID']] = scores.get(r['StructureID'], [None, None])
                scores[r['StructureID']][0] = -r['total']
                wts += r['total']
            elif r['ScoreType'] == 'Mutant':
                scores[r['StructureID']] = scores.get(r['StructureID'], [None, None])
                scores[r['StructureID']][1] = r['total']
                muts += r['total']
        for k, v in scores.iteritems():
            print(k, v)
            scores[k] = sum(v)
        for k, v in sorted(scores.iteritems(), key = lambda x:x[1]):
            print(k, v)
        print(len(scores))
        print(wts/len(scores))
        print(muts/len(scores))


        sys.exit(0)
        dumpfile = ddG_connection.getData(72856)
        print(len(dumpfile))
        write_file('72856.zip', dumpfile, 'wb')
        sys.exit(0)

    x = True

    prediction_set_for_analysis = 'Protocol_16_r57471'




    #for analysis_subset in ['Kellogg', 'Guerois', 'Potapov', 'CuratedProTherm']:
    for analysis_subset in ['CuratedProTherm']:
        for prediction_set_for_analysis in ['RosCon2013_P16_score12prime', 'RosCon2013_P16_talaris2013', 'RosCon2013_P16_talaris2013sc', 'Protocol_16_r57471']:
            if prediction_set_for_analysis == 'RosCon2013_P16_score12prime':
                analysis_breakdown = load_data('score12_analysis_input.json')
            elif prediction_set_for_analysis == 'RosCon2013_P16_talaris2013':
                analysis_breakdown = load_data('r55534_talaris2013_analysis_input.json')
            elif prediction_set_for_analysis == 'RosCon2013_P16_talaris2013sc':
                analysis_breakdown = load_data('r55534_talaris2013sc_analysis_input.json')
            elif prediction_set_for_analysis == 'Protocol_16_r57471':
                analysis_breakdown = load_data('r57471_analysis_input.json')

            colortext.message("\n\n\n\n*** %s: %s ***" % (analysis_subset, prediction_set_for_analysis))
        #for analysis_subset in analysis_breakdown.analysis_datasets.keys():

            #if analysis_subset != 'Kellogg':
            #    continue
            if analysis_subset not in ['Kellogg', 'AlaScan-GPK', 'TransmembraneProteins', 'CuratedProTherm', 'Guerois', 'Potapov']:
                continue
            #if analysis_subset not in ['Kellogg']:
            #    continue
            if analysis_subset not in ['Kellogg', 'Guerois', 'Potapov', 'CuratedProTherm', ]:
                continue
            ddGdb = ddgdbapi.ddGDatabase()
            rrecords = ddGdb.execute_select('''
    SELECT DISTINCT Prediction.ID AS PredictionID, UserAnalysisSet.RecordNumber
    FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperimentID=UserDataSetExperiment.ID
    INNER JOIN UserAnalysisSet ON UserDataSetExperiment.ExperimentID=UserAnalysisSet.ExperimentID
    WHERE PredictionSet=%s
    AND UserDataSetExperiment.UserDataSetID=UserAnalysisSet.UserDataSetID
    and UserDataSetExperiment.PDBFileID=UserAnalysisSet.PDB_ID
    AND Subset=%s
    ORDER BY  RecordNumber''', parameters=(prediction_set_for_analysis, analysis_subset,))
            pred_record_map = {}
            for rrecord in rrecords:
                pred_record_map[rrecord['PredictionID']] = [rrecord['RecordNumber']]

            #for scoring_method in ['Kellogg', 'Noah']:
            for scoring_method in ['Kellogg_top3']:
                print('Analyzing %d records.' % len(rrecords))
                analysis_breakdown.analyze_subset_all(analysis_subset, scoring_method, pred_record_map)
                analysis_breakdown.analyze_subset_by_mutation_size(analysis_subset, scoring_method)
                continue
                analysis_breakdown.analyze_subset_single_no_GP(analysis_subset, scoring_method)
                analysis_breakdown.analyze_subset_single_GP(analysis_subset, scoring_method)
                analysis_breakdown.analyze_subset_multiple(analysis_subset, scoring_method)

                #analysis_breakdown.analyze_subset_by_secondary_structure(analysis_subset, scoring_method)
                #analysis_breakdown.analyze_subset_by_binned_chain_length(analysis_subset, scoring_method, num_bins = 9)

                #analysis_breakdown.analyze_subset_by_aromaticity(analysis_subset, scoring_method)
                #analysis_breakdown.analyze_subset_by_exposure(analysis_subset, scoring_method, cut_off=0.4)
                #analysis_breakdown.analyze_subset_by_wildtype_charge(analysis_subset, scoring_method)
                print('')
                #analysis_breakdown.analyze_subset_by_specific_resolutions(analysis_subset, scoring_method)
                #analysis_breakdown.analyze_subset_by_binned_resolutions(analysis_subset, scoring_method)
            #break
