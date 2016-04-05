prediction_set_id = 'james-backrub-rscript-full'

user_dataset_name = 'AllBindingAffinity'
tagged_subset = 'ZEMu'
prediction_set_description = "James' first attempt to run backrub ddG with a Rosetta script, trials set to 1000, nstruct 10, kt=1.6"
extra_flags = []

score_method_ids = [12]

# Prediction set analysis
prediction_set_credit = 'James Lucas'
take_lowests = [3, 10]
expectn = 10
use_existing_benchmark_data = False
allow_missing_case_failures = True

temporary_data_location = '/kortemmelab/home/james.lucas/160322-backrub-full-reduce-MC-Trials/output'
