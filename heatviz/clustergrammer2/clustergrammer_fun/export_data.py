# def export_net_json(net, net_type, indent='no-indent'):
#   ''' export json string of dat '''
#   import json
#   from copy import deepcopy

#   if net_type == 'dat':
#     exp_dict = deepcopy(net.dat)

#     if type(exp_dict['mat']) is not list:
#       exp_dict['mat'] = exp_dict['mat'].tolist()

#   elif net_type == 'viz':
#     exp_dict = net.viz

#   elif net_type == 'sim_row':
#     exp_dict = net.sim['row']

#   elif net_type == 'sim_col':
#     exp_dict = net.sim['col']

#   # make json
#   if indent == 'indent':
#     exp_json = json.dumps(exp_dict, indent=2)
#   else:
#     exp_json = json.dumps(exp_dict)

#   return exp_json

import json

class RoundingEncoder(json.JSONEncoder):
    """Custom JSON encoder that rounds floats to 3 decimal places"""
    def default(self, obj):
        if isinstance(obj, float):
            return round(obj, 3)
        return super().default(obj)

    def encode(self, obj):
        # Override encode to handle floats in nested structures
        return super().encode(self._round_floats(obj))

    def _round_floats(self, obj):
        if isinstance(obj, float):
            return round(obj, 3)
        elif isinstance(obj, dict):
            return {k: self._round_floats(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._round_floats(item) for item in obj]
        return obj


def export_net_json(net, net_type, indent='no-indent'):
    ''' export json string of dat '''
    import json
    import numpy as np
    from copy import deepcopy

    if net_type == 'dat':
        exp_dict = deepcopy(net.dat)

        if type(exp_dict['mat']) is not list:
            exp_dict['mat'] = exp_dict['mat'].tolist()

    elif net_type == 'viz':
        exp_dict = net.viz

    elif net_type == 'sim_row':
        exp_dict = net.sim['row']

    elif net_type == 'sim_col':
        exp_dict = net.sim['col']

    # Clean NaN values
    exp_dict = clean_nan_for_json(exp_dict)

    # make json with custom encoder that rounds floats
    if indent == 'indent':
        exp_json = json.dumps(exp_dict, indent=2, cls=RoundingEncoder)
    else:
        exp_json = json.dumps(exp_dict, cls=RoundingEncoder)

    return exp_json

def clean_nan_for_json(obj):
    """
    Recursively clean NaN, inf, and -inf values from nested data structures.
    Converts them to None (which becomes null in JSON).
    """
    import numpy as np
    
    if isinstance(obj, dict):
        return {k: clean_nan_for_json(v) for k, v in obj.items()}
    
    elif isinstance(obj, list):
        return [clean_nan_for_json(item) for item in obj]
    
    elif isinstance(obj, np.ndarray):
        # Convert numpy array to list and clean
        return clean_nan_for_json(obj.tolist())
    
    elif isinstance(obj, (float, np.floating)):
        if np.isnan(obj):
            return None  # NaN → null in JSON
        elif np.isinf(obj):
            if obj > 0:
                return "Infinity"
            else:
                return "-Infinity"
        else:
            return float(obj)  # Convert numpy float to Python float
    
    elif isinstance(obj, (int, np.integer)):
        return int(obj)  # Ensure it's a Python int, not numpy int
    
    elif isinstance(obj, (np.bool_, bool)):
        return bool(obj)  # Ensure it's a Python bool
    
    else:
        return obj


def write_matrix_to_tsv(net, filename=None, df=None):
  '''
  This will export the matrix in net.dat or a dataframe (optional df in
  arguments) as a tsv file. Row/column categories will be saved as tuples in
  tsv, which can be read back into the network object.
  '''
  import pandas as pd

  if df is None:
    df = net.dat_to_df()

  return df.to_csv(filename, sep='\t')

def write_json_to_file(net, net_type, filename, indent='no-indent'):

  exp_json = net.export_net_json(net_type, indent)

  fw = open(filename, 'w')
  fw.write(exp_json)
  fw.close()

def save_dict_to_json(inst_dict, filename, indent='no-indent'):
  import json
  fw = open(filename, 'w')
  if indent == 'indent':
    fw.write(json.dumps(inst_dict, indent=2))
  else:
    fw.write(json.dumps(inst_dict))
  fw.close()