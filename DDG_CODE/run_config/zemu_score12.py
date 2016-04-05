prediction_set_id = 'ddgmon16-zemu-score12' # This should be a valid filename
protocol_name = 'ddgmon16-score12'

user_dataset_name = 'AllBindingAffinity'
tagged_subset = 'ZEMu lite'
prediction_set_description = 'New run of ddg_monomer_16 on ZEMU, hopefully with occupancy/residue loading issues resolved'
extra_flags = ['-restore_pre_talaris2013_behavior', '-score:patch score12.wts_patch']

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
