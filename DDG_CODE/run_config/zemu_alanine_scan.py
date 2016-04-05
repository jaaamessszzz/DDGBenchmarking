prediction_set_id = 'alascan-zemu-talaris2014' # This should be a valid filename
protocol_name = 'alascan'

rosetta_scripts_command_line = '-ignore_zero_occupancy false -ignore_unrecognized_res'
rosetta_scripts_file = 'interface/alanine_scanning.xml'
rosetta_scripts_score_fxn = 'interface'

user_dataset_name = 'AllBindingAffinity'
tagged_subset = 'ZEMu'
prediction_set_description = 'alanine scanning on ZEMU with default (Talaris 2014) score function'
extra_flags = ['-ignore_zero_occupancy false', '-ignore_unrecognized_res']

# score_method_id = 7 # rescore with interface weights score method (this is probably best...need to benchmark more)
# score_method_id = 8 # rescore with talaris weights score method
# score_method_ids = [7, 8, 9, 10]

# Prediction set analysis
prediction_set_credit = 'Kyle Barlow'
take_lowests = [1]
expectn = 1
