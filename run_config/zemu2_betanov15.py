prediction_set_id = 'ddg_monomer_16-zemu-betanov15' # This should be a valid filename
protocol_name = 'ddg_monomer_16-betanov15'

user_dataset_name = 'AllBindingAffinity'
tagged_subset = 'ZEMu'
prediction_set_description = 'ddg_monomer_16 on ZEMU with beta Nov 2015 score function'
extra_flags = ['-beta_nov15']

# score_method_id = 7 # rescore with interface weights score method (this is probably best...need to benchmark more)
# score_method_id = 8 # rescore with talaris weights score method
score_method_ids = [10]

# Prediction set analysis
prediction_set_credit = 'Kyle Barlow'
take_lowests = [3, 10, 20, 30, 40, 50]
expectn = 50
