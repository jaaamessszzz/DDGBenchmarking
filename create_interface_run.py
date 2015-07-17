import os, sys

# Add parent directory to path
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from ddglib.ppi_api import get_interface_with_config_file
    
prediction_set_id = 'ZEMu run 1'
ppi_api = get_interface_with_config_file()

ppi_api.help()

# Create the prediction set
#prediction_set_id = 'ZEMu run 1 with my protocol, r57982'
#ppi_api.add_prediction_set(prediction_set_id, halted = True, priority = 7, batch_size = 41, allow_existing_prediction_set = True)

# Unnecessary but here is how to change the values of batch_size, priority
# ppi_api.alter_prediction_set_batch_size(prediction_set_id, 40)
# ppi_api.alter_prediction_set_priority(prediction_set_id, 5)

# This should be called before kicking off jobs (or set halted = False above)
#ppi_api.start_prediction_set(prediction_set_id)

# compile the python submission script
