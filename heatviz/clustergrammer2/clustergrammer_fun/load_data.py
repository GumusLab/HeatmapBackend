import io, sys
import json
import pandas as pd
import numpy as np
from . import categories
from . import proc_df_labels
from . import data_formats
from . import make_unique_labels

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

def load_file(net, filename):
  # reset network when loaing file, prevents errors when loading new file
  # have persistent categories

  net.reset()

  f = open(filename, 'r')

  file_string = f.read()
  f.close()

  load_file_as_string(net, file_string, filename)

def load_file_as_string(net, file_string, filename=''):

  if (sys.version_info > (3, 0)):
    # python 3
    ####################
    file_string = str(file_string)
  else:
    # python 2
    ####################
    file_string = unicode(file_string)

  buff = io.StringIO(file_string)

  if '/' in filename:
    filename = filename.split('/')[-1]

  net.load_tsv_to_net(buff, filename)

def load_stdin(net):
  data = ''

  for line in sys.stdin:
    data = data + line

  data = StringIO.StringIO(data)

  net.load_tsv_to_net(data)

# def load_tsv_to_net(net, file_buffer, filename=None):
#   lines = file_buffer.getvalue().split('\n')
#   num_labels = categories.check_categories(lines)

#   print(' ******* num labels are as follows *****', num_labels)

#   row_arr = list(range(num_labels['row']))
#   col_arr = list(range(num_labels['col']))

#   print(row_arr)
#   print(col_arr)

#   # use header if there are col categories
#   if len(col_arr) > 1:
#     df = pd.read_table(file_buffer, index_col=row_arr,
#                                   header=col_arr)
#   else:
#     df = pd.read_table(file_buffer, index_col=row_arr)

#   df = proc_df_labels.main(df)

#   net.df_to_dat(df, True)
#   net.dat['filename'] = filename

def load_tsv_to_net(net, file_buffer, filename=None):
    lines = file_buffer.getvalue().split('\n')
    num_labels = categories.check_categories(lines)

    print(' ******* num labels are as follows *****', num_labels)

    row_arr = list(range(num_labels['row']))
    col_arr = list(range(num_labels['col']))

    print(row_arr)
    print(col_arr)

    # use header if there are col categories
    if len(col_arr) > 1:
        df = pd.read_table(file_buffer, index_col=row_arr, header=col_arr)
    else:
        df = pd.read_table(file_buffer, index_col=row_arr)

    # 🔧 FIX: Convert all data to numeric, handle mixed types
    print("🧹 Cleaning data types...")
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Remove rows/columns with all NaN
    df = df.dropna(how='all').dropna(axis=1, how='all')
    
    # 🔧 CRITICAL FIX: Ensure index and columns are strings
    print("🔧 Converting index and columns to strings...")
    if isinstance(df.index, pd.MultiIndex):
      df.index = df.index.astype(object)
    else:
        df.index = df.index.astype(str)

    if isinstance(df.columns, pd.MultiIndex):
        df.columns = df.columns.astype(object) 
    else:
        df.columns = df.columns.astype(str)
    
    print(f"✅ Final data shape: {df.shape}")
    print(f"✅ Index type: {type(df.index[0])}")
    print(f"✅ Column type: {type(df.columns[0])}")

    df = proc_df_labels.main(df)
    net.df_to_dat(df, True)
    net.dat['filename'] = filename


# def load_tsv_to_net(net, file_buffer, filename=None):
#     """Fixed version that handles MultiIndex columns properly."""
    
#     lines = file_buffer.getvalue().split('\n')
#     num_labels = categories.check_categories(lines)
    
#     print(f"🔍 TSV LOADING DEBUG:")
#     print(f"   Total lines: {len(lines)}")
#     print(f"   Category detection: {num_labels}")
    
#     row_arr = list(range(num_labels['row']))
#     col_arr = list(range(num_labels['col']))
    
#     print(f"   Row headers (index_col): {row_arr}")
#     print(f"   Col headers (header): {col_arr}")
    
#     # Reset buffer before pandas reads it
#     file_buffer.seek(0)
    
#     print(f"🔍 Pandas loading...")
#     try:
#         if len(col_arr) > 1:
#             df = pd.read_table(file_buffer, index_col=row_arr, header=col_arr)
#         else:
#             df = pd.read_table(file_buffer, index_col=row_arr)
            
#         print(f"   DataFrame shape after pandas: {df.shape}")
#         print(f"   Column type: {type(df.columns)}")
        
#         # 🔧 CRITICAL FIX: Handle MultiIndex columns
#         if isinstance(df.columns, pd.MultiIndex):
#             print(f"   ⚠️ MultiIndex columns detected with {len(col_arr)} levels")
#             print(f"   MultiIndex levels: {df.columns.nlevels}")
#             print(f"   Sample columns: {df.columns[:3].tolist()}")
            
#             # Flatten the MultiIndex columns
#             print(f"   🔧 Flattening MultiIndex columns...")
            
#             # Create flattened column names by joining all levels
#             flattened_columns = []
#             for col_tuple in df.columns:
#                 # Join non-empty parts of the tuple with underscore
#                 col_parts = [str(part) for part in col_tuple if str(part) != 'nan' and str(part).strip() != '']
#                 if col_parts:
#                     flattened_name = '_'.join(col_parts)
#                 else:
#                     # If all parts are empty/nan, use a default name
#                     flattened_name = f'col_{len(flattened_columns)}'
#                 flattened_columns.append(flattened_name)
            
#             # Assign flattened column names
#             df.columns = flattened_columns
#             print(f"   ✅ Columns flattened: {len(flattened_columns)} columns")
#             print(f"   Sample flattened columns: {df.columns[:5].tolist()}")
        
#         # Ensure column names are strings and unique
#         df.columns = df.columns.astype(str)
        
#         # Handle duplicate column names
#         if df.columns.duplicated().any():
#             print(f"   ⚠️ Duplicate column names found, making unique...")
#             df.columns = pd.Index([f"{col}_{i}" if df.columns[:i+1].duplicated()[i] 
#                                  else col for i, col in enumerate(df.columns)])
#             print(f"   ✅ Column names made unique")
        
#     except Exception as e:
#         print(f"❌ Pandas loading failed: {e}")
#         import traceback
#         traceback.print_exc()
#         return
    
#     print(f"🔍 Label processing...")
#     df_before_proc = df.copy()
#     df = proc_df_labels.main(df)
#     print(f"   Shape before proc_df_labels: {df_before_proc.shape}")
#     print(f"   Shape after proc_df_labels: {df.shape}")
    
#     print(f"🔍 Network conversion...")
#     try:
#         # Initialize network structure
#         if not hasattr(net, 'dat') or net.dat is None:
#             net.dat = {}
#         if 'nodes' not in net.dat:
#             net.dat['nodes'] = {}
#         if 'node_info' not in net.dat:
#             net.dat['node_info'] = {'row': {}, 'col': {}}
#         if not hasattr(net, 'meta_cat'):
#             net.meta_cat = False
        
#         print(f"   Pre-conversion DataFrame info:")
#         print(f"     Shape: {df.shape}")
#         print(f"     Index type: {type(df.index)}")
#         print(f"     Columns type: {type(df.columns)}")
#         print(f"     Index sample: {df.index[:3].tolist()}")
#         print(f"     Columns sample: {df.columns[:3].tolist()}")
        
#         # Call the df_to_dat function
#         net.df_to_dat(df, True)
#         net.dat['filename'] = filename
        
#         print(f"   ✅ Network conversion successful")
        
#         # Verify the result
#         if 'mat' in net.dat:
#             final_shape = np.array(net.dat['mat']).shape
#             print(f"   Final network matrix shape: {final_shape}")
#             print(f"   node_info row length: {len(net.dat['node_info']['row'])}")
#             print(f"   node_info col length: {len(net.dat['node_info']['col'])}")
        
#     except Exception as e:
#         print(f"❌ Network conversion failed: {e}")
#         import traceback
#         traceback.print_exc()
#         raise


def debug_categories_check_categories(lines):
    """Debug the categories.check_categories function."""
    print(f"🔍 DEBUGGING categories.check_categories:")
    print(f"   Input: {len(lines)} lines")
    print(f"   First 3 lines:")
    for i, line in enumerate(lines[:3]):
        print(f"      Line {i}: {repr(line[:100])}")
    
    # Call the actual function
    result = categories.check_categories(lines)
    print(f"   Result: {result}")
    
    return result







def load_json_to_dict(filename):
  f = open(filename, 'r')
  inst_dict = json.load(f)
  f.close()
  return inst_dict

def load_gmt(filename):
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  gmt = {}
  for i in range(len(lines)):
    inst_line = lines[i].rstrip()
    inst_term = inst_line.split('\t')[0]
    inst_elems = inst_line.split('\t')[2:]
    gmt[inst_term] = inst_elems

  return gmt

def load_data_to_net(net, inst_net):
  ''' load data into nodes and mat, also convert mat to numpy array'''
  net.dat['nodes'] = inst_net['nodes']
  net.dat['mat'] = inst_net['mat']
  data_formats.mat_to_numpy_arr(net)