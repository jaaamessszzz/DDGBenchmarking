prediction_set_id = 'ran_gsp1_ddg-monomer-16_1' # This should be a valid filename
protocol_name = 'ddg_monomer_16_003'

user_dataset_name = 'RAN-GSP'
tagged_subset = None
prediction_set_description = 'First run of ddg monomer (row 16) on ran gsp1 data'
extra_flags = ['-in:auto_setup_metals']
keep_hetatm_lines = True
keep_all_lines = True

# score_method_id = 7 # rescore with interface weights score method (this is better for this run)
# score_method_id = 8 # rescore with talaris weights score method (not as good for this run)
score_method_ids = [7, 8]

# Prediction set analysis
prediction_set_credit = 'Kyle Barlow'
# take_lowests = range(1,51)
take_lowests = [3, 50]
expectn = 50
use_existing_benchmark_data = True
allow_missing_case_failures = False
