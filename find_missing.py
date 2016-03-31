import sys, os

from ddglib.ppi_api import get_interface_with_config_file
import klab.cluster_template.parse_settings as parse_settings
from ddglib.ddg_monomer_ppi_api import get_interface as get_interface_factory
from create_ddg_monomer_interface_run import main as create_ddg_monomer_interface_run
import imp

def main():
    assert( len(sys.argv) >= 1 )

    cfg = imp.load_source('cfg', sys.argv[1])

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']
    ppi_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = settings['local_rosetta_installation_path'] + '/database', get_interface_factory = get_interface_factory )

    prediction_set_id = cfg.prediction_set_id

    prediction_ids = sorted(ppi_api.get_prediction_ids(prediction_set_id))

    results = ppi_api.DDG_db.execute_select('''
    SELECT
    PredictionPPIStructureScore.PredictionPPIID, PredictionPPIStructureScore.ScoreMethodID,
    PredictionPPIStructureScore.StructureID,
    COUNT(PredictionPPIStructureScore.ScoreType) AS ScoreTypeCount
    FROM PredictionPPIStructureScore
    INNER JOIN PredictionPPI ON PredictionPPI.ID=PredictionPPIStructureScore.PredictionPPIID
    WHERE PredictionPPI.PredictionSet='%s'
    GROUP BY
    PredictionPPIStructureScore.PredictionPPIID,
    PredictionPPIStructureScore.ScoreMethodID,
    PredictionPPIStructureScore.StructureID
    ''' % prediction_set_id)
    structure_ids = set()
    count_dict = {}
    expected_structures_counts = {}
    score_methods = set()
    for row in results:
        id_tup = (row['PredictionPPIID'], row['ScoreMethodID'])
        if id_tup not in count_dict:
            count_dict[id_tup] = {}
        count_dict[id_tup][row['StructureID']] = row['ScoreTypeCount']
        structure_ids.add(row['StructureID'])
        if row['ScoreTypeCount'] not in expected_structures_counts:
            expected_structures_counts[ row['ScoreTypeCount'] ] = 0
        expected_structures_counts[ row['ScoreTypeCount'] ] += 1
        score_methods.add( row['ScoreMethodID'] )
    expected_structure_ids = set(range(1, max(structure_ids) + 1))
    missing_structure_ids = expected_structure_ids.difference(structure_ids)
    if len(missing_structure_ids) > 0:
        print 'Missing structure ids from all records:', sorted(missing_structure_ids)
    expected_structures = sorted([(v,k) for k,v in expected_structures_counts.iteritems()])[-1][1]
    ids_to_repeat = set()
    for score_method in sorted(score_methods):
        for prediction_id in prediction_ids:
            id_tup = prediction_id, score_method
            if id_tup not in count_dict:
                print id_tup, 'missing all scores'
                ids_to_repeat.add( prediction_id )
            else:
                for struct_id in sorted(structure_ids):
                    if struct_id not in count_dict[id_tup]:
                        print id_tup, 'missing scores for structure', struct_id
                        ids_to_repeat.add( prediction_id )
                    elif count_dict[id_tup][struct_id] != expected_structures:
                        print id_tup, 'missing some scores for structure', struct_id, 'has only', count_dict[id_tup][struct_id], '(expecting: ', expected_structures, ')'
                        ids_to_repeat.add( prediction_id )
    print 'IDs to repeat:', ids_to_repeat
    if len( ids_to_repeat ) > 0:
        create_ddg_monomer_interface_run( prediction_ids = sorted(ids_to_repeat), memory_free = '6.0G', cfg = cfg )
if __name__ == '__main__':
    main()
