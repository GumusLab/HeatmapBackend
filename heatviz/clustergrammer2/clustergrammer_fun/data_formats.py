from . import make_unique_labels
import pandas as pd
from . import categories

def df_to_dat(net, df, define_cat_colors=False):
  '''
  This is always run when data is loaded.
  '''

  # check if df has unique values
  df = make_unique_labels.main(net, df)

  net.dat['mat'] = df.values
  net.dat['nodes']['row'] = df.index.tolist()
  net.dat['nodes']['col'] = df.columns.tolist()

  # if net.meta_cat == False or net.is_downsampled:
  if net.meta_cat == False:

    # tuple cats
    ##################################

    for axis in ['row', 'col']:
      
      inst_nodes = net.dat['nodes'][axis]

      # if type(inst_nodes[0]) is tuple:

      #    # 🔍 DEBUG CODE START
      #   print(f"🔍 TUPLE DEBUGGING for {axis}:")
      #   print(f"   First tuple: {inst_nodes[0]}")
      #   print(f"   Tuple length: {len(inst_nodes[0])}")
      #   num_cats = len(inst_nodes[0]) - 1
      #   print(f"   Expected categories: {num_cats}")
      #   print(f"   Will try to access indices: {list(range(1, num_cats + 1))}")
        
      #   lengths = [len(x) for x in inst_nodes[:5]]
      #   print(f"   First 5 tuple lengths: {lengths}")
        
      #   if len(set(lengths)) > 1:
      #       print(f"   ❌ ERROR: Tuples have different lengths!")
      #   # 🔍 DEBUG CODE END

      #   if axis == 'row':
      #     net.dat['node_info'][axis]['full_names'] = df.index.tolist()
      #   elif axis == 'col':
      #     net.dat['node_info'][axis]['full_names'] = df.columns.tolist()

      #   # get the number of categories from the length of the tuple
      #   # subtract 1 because the name is the first element of the tuple
      #   num_cats = len(inst_nodes[0]) - 1
      #   for inst_cat in range(num_cats):
      #     cat_name = 'cat-' + str(inst_cat)
      #     cat_index = inst_cat + 1
      #     cat_values = [x[cat_index] for x in inst_nodes]
      #     net.dat['node_info'][axis][cat_name] = cat_values

      #   # clean up nodes after parsing categories
      #   net.dat['nodes'][axis] = [x[0] for x in inst_nodes]
      

      
      if type(inst_nodes[0]) is tuple:
    
          num_cats = len(inst_nodes[0]) - 1
      
          
          lengths = [len(x) for x in inst_nodes[:5]]
          
          if len(set(lengths)) > 1:
              print(f"ERROR: Tuples have different lengths!")
          
          # 🔍 ADD MORE DEBUGGING HERE
          try:
              if axis == 'row':
                  net.dat['node_info'][axis]['full_names'] = df.index.tolist()
              elif axis == 'col':
                  net.dat['node_info'][axis]['full_names'] = df.columns.tolist()
              
              for inst_cat in range(num_cats):
                  cat_name = 'cat-' + str(inst_cat)
                  cat_index = inst_cat + 1
                  
                  try:
                      cat_values = [x[cat_index] for x in inst_nodes]
                      
                      # Check if all values are NaN
                      nan_count = sum(1 for val in cat_values if pd.isna(val) or str(val) == 'nan')
                      total_count = len(cat_values)
                      
                      if nan_count == total_count:
                          continue  # Skip this category entirely
                      else:
                          net.dat['node_info'][axis][cat_name] = cat_values
                  except Exception as cat_e:
                      print(f"   ❌ ERROR in category {cat_name}: {cat_e}")
                      raise
              
              print(f"   🔧 Cleaning up nodes...")
              try:
                  net.dat['nodes'][axis] = [x[0] for x in inst_nodes]
                  print(f"   ✅ Nodes cleaned up: {len(net.dat['nodes'][axis])} items")
              except Exception as cleanup_e:
                  print(f"ERROR in cleanup: {cleanup_e}")
                  raise
                  
          except Exception as e:
              print(f" ERROR in tuple processing for {axis}: {e}")
              import traceback
              traceback.print_exc()
              raise
      
      else:
          # Handle non-tuple nodes (simple strings)
          
          if axis == 'row':
              net.dat['node_info'][axis]['full_names'] = df.index.tolist()
          elif axis == 'col':
              net.dat['node_info'][axis]['full_names'] = df.columns.tolist()          
    # Continue with your existing code after the if block...

  else:
    

    for axis in ['row', 'col']:

      inst_nodes = net.dat['nodes'][axis]

      if axis == 'row':
        net.dat['node_info'][axis]['full_names'] = df.index.tolist()
      elif axis == 'col':
        net.dat['node_info'][axis]['full_names'] = df.columns.tolist()

      inst_cats = []
      if axis == 'row':
        # inst_cats = net.meta_row.columns.tolist()
        if hasattr(net, 'row_cats'):
          inst_cats = net.row_cats
      else:
        # inst_cats = net.meta_col.columns.tolist()
        if hasattr(net, 'col_cats'):
          inst_cats = net.col_cats

      num_cats = len(inst_cats)

      # if axis == 'row':
      #   num_cats = len(net.row_cats)
      # elif axis == 'col':
      #   num_cats = len(net.col_cats)

      if num_cats > 0:
        for inst_cat in range(num_cats):
          cat_name = 'cat-' + str(inst_cat)
          cat_index = inst_cat + 1

          cat_title = inst_cats[inst_cat]
          if axis == 'row':
            if net.is_downsampled:
              if hasattr(net, 'meta_ds_row'):
                cat_values = net.meta_ds_row.loc[inst_nodes, cat_title].apply(lambda x: cat_title + ': ' + str(x)).values.tolist()
              else:
                cat_values = net.meta_row.loc[inst_nodes, cat_title].apply(lambda x: cat_title + ': ' + str(x)).values.tolist()

            # detault with no downsampling
            else:
              cat_values = net.meta_row.loc[inst_nodes, cat_title].apply(lambda x: cat_title + ': ' + str(x)).values.tolist()
          else:
            # cat_values = net.meta_col.loc[inst_nodes, cat_title].apply(lambda x: cat_title + ': ' + x).values.tolist()
            if net.is_downsampled:
              if hasattr(net, 'meta_ds_col'):
                # print(inst_nodes)

                cat_values = net.meta_ds_col.loc[inst_nodes, cat_title].apply(lambda x: cat_title + ': ' + str(x)).values.tolist()
              else:
                cat_values = net.meta_col.loc[inst_nodes, cat_title].apply(lambda x: cat_title + ': ' + str(x)).values.tolist()

            # detault with no downsampling
            else:
              cat_values = net.meta_col.loc[inst_nodes, cat_title].apply(lambda x: cat_title + ': ' + str(x)).values.tolist()

          net.dat['node_info'][axis][cat_name] = cat_values


  for axis in ['row', 'col']:
      if axis in net.dat.get('node_info', {}):
          node_info = net.dat['node_info'][axis]
          print(f"     {axis}: {list(node_info.keys())}")
          for key, values in node_info.items():
              if isinstance(values, list):
                  print(f"       {key}: {len(values)} items")
              else:
                  print(f"       {key}: {type(values)}")

  try:
      categories.dict_cat(net, define_cat_colors=define_cat_colors)
  except Exception as e:
     
      import traceback
      traceback.print_exc()
      raise

  # categories.dict_cat(net, define_cat_colors=define_cat_colors)

def dat_to_df(net):

  nodes = {}
  for axis in ['row', 'col']:
    if 'full_names' in net.dat['node_info'][axis]:
      nodes[axis] = net.dat['node_info'][axis]['full_names']
    else:
      nodes[axis] = net.dat['nodes'][axis]

  df = pd.DataFrame(data=net.dat['mat'], columns=nodes['col'],
      index=nodes['row'])

  return df

def mat_to_numpy_arr(self):
  ''' convert list to numpy array - numpy arrays can not be saved as json '''
  import numpy as np
  self.dat['mat'] = np.asarray(self.dat['mat'])