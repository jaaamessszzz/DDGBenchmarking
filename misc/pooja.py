import sys
import pprint

if False:
  if __name__ == "__main__":
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, "../updatedb")
    sys.path.insert(0, '/home/oconchus/dev/')
    sys.path.insert(0, "/home/oconchus/dev/klab")

import klab.colortext as colortext
from ddglib.ppi_api import get_interface_with_config_file as get_ppi_interface_with_config_file


# Set up database connection
try:
    ppi_api = get_ppi_interface_with_config_file(host_config_name = 'kortemmelab')
except:
    colortext.error('Database connection failed.')
    raise
colortext.message('Connected to database.')


# Pick a scoring method
score_method_id = ppi_api.get_score_method_id('Rescore-Talaris2014', method_authors = 'kyle', method_type = 'ddg_monomer rescore')

# Get the best structures for prediction 23849
pprint.pprint(ppi_api.get_top_x_scores(23849, score_method_id, 'WildTypeComplex', 3, component = 'total', order_by = 'ASC'))
pprint.pprint(ppi_api.get_top_x_scores(23849, score_method_id, 'MutantComplex', 3, component = 'total', order_by = 'ASC'))

