import os, sys

import klab.cluster_template.parse_settings as parse_settings
from ddglib.ppi_api import get_interface_with_config_file
from ddglib.api_api import get_interface as get_interface_factory
from jameslucas.UtilityScripts.RosettaOut_Parser import parse_rosetta_out
import importlib

def main():
    assert( len(sys.argv) >= 2 )
    import_module = sys.argv[1]
    if import_module.endswith('.py'):
        import_module = import_module[:-3]
    if '/' in import_module:
        import_module = import_module.replace('/', '.')
    cfg = importlib.import_module(import_module, package=None)

    prediction_set_name = cfg.prediction_set_id

    settings = parse_settings.get_dict()
    rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

    api_api = get_interface_with_config_file(rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = '/home/kyleb/rosetta/working_branches/alascan/database', get_interface_factory = get_interface_factory )
    score_method_ids = cfg.score_method_ids

    scores_dict, struct_count_dict = parse_rosetta_out(cfg.temporary_data_location, verbose = False)
    api_api.set_run_data( scores_dict )

    for score_method_id in score_method_ids:
        if len(score_method_ids) > 1:
            print 'Processing score_method_id: %d\n' % score_method_id
        api_api.extract_data(prediction_set_name, score_method_id = score_method_id, expectn = cfg.expectn, root_directory = cfg.temporary_data_location )


if __name__ == '__main__':
    main()
