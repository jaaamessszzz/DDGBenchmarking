import sys
import pprint
import glob

if True:
  if __name__ == "__main__":
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, "../updatedb")
    sys.path.insert(0, '/home/oconchus/dev/')
    sys.path.insert(0, "/home/oconchus/dev/klab")
else:
    import klab

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
wild_type_complexes = ppi_api.get_top_x_scores(23849, score_method_id, 'WildTypeComplex', 3, component = 'total', order_by = 'ASC')
wild_type_filenames = []
for wtc in wild_type_complexes:
    wild_type_filenames.append([f for f in glob.glob('repacked_wt*_round_{0}.*'.format(wtc['StructureID']))][0])
print(wild_type_filenames)


mutant_complexes = ppi_api.get_top_x_scores(23849, score_method_id, 'MutantComplex', 3, component = 'total', order_by = 'ASC')
mutant_filenames = []
for mtc in mutant_complexes:
    mutant_filenames.append([f for f in glob.glob('mut*_round_{0}.*'.format(wtc['StructureID']))][0])
print(mutant_filenames)

