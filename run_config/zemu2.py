prediction_set_id = 'ddg_monomer_16_003-zemu-2' # This should be a valid filename
protocol_name = 'ddg_monomer_16_003'

user_dataset_name = 'AllBindingAffinity'
tagged_subset = 'ZEMu'
prediction_set_description = 'New run of ddg_monomer_16 on ZEMU, hopefully with occupancy/residue loading issues resolved'
extra_flags = []

# score_method_id = 7 # rescore with interface weights score method (this is better for this run)
# score_method_id = 8 # rescore with talaris weights score method (not as good for this run)
score_method_ids = [7, 8]

# Prediction set analysis
prediction_set_credit = 'Kyle Barlow'
take_lowests = range(1,51)
expectn = 50
