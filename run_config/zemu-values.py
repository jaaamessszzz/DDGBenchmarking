prediction_set_id = 'zemu-values' # This should be a valid filename
protocol_name = 'zemu-paper'

user_dataset_name = 'AllBindingAffinity'
tagged_subset = 'ZEMu'
prediction_set_description = 'ZEMu author values'
extra_flags = []

# score_method_id = 7 # rescore with interface weights score method (this is better for this run)
# score_method_id = 8 # rescore with talaris weights score method (not as good for this run)
score_method_ids = [11]

# Prediction set analysis
prediction_set_credit = 'ZEMu authors'
# take_lowests = range(1,51)
take_lowests = [1]
expectn = 1
use_existing_benchmark_data = True
allow_missing_case_failures = False
