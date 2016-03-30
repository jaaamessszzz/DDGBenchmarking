backrub_temp = 2.4
inner_jobs = 50 # Number of structures
prediction_set_id = 'zemu-psbrub_%.1f' % backrub_temp # This should be a valid filename

user_dataset_name = 'AllBindingAffinity'
tagged_subset = 'ZEMu lite'
prediction_set_description = 'Backrub ddG (Pooja and Samuel) on ZEMu with kT of %.1f' % backrub_temp
extra_flags = []

# score_method_id = 7 # rescore with interface weights score method (this is probably best...need to benchmark more)
# score_method_id = 8 # rescore with talaris weights score method
score_method_ids = [8]

# Prediction set analysis
prediction_set_credit = 'Kyle Barlow'
take_lowests = [3, 50]#, 10, 20, 30, 40, 50]
expectn = 40
use_existing_benchmark_data = False
allow_missing_case_failures = True
