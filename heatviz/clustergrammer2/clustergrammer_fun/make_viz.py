def viz_json(net, dendro=True, links=False):
  ''' make the dictionary for the clustergram.js visualization '''
  from . import calc_clust
  import numpy as np

  # linkage information
  # net.viz['linkage'] = {}
  # net.viz['linkage']['row'] = net.dat['node_info']['row']['Y'].tolist()
  # net.viz['linkage']['col'] = net.dat['node_info']['col']['Y'].tolist()

  all_dist = calc_clust.group_cutoffs()


  # node information
  for inst_rc in net.dat['nodes']:

    inst_keys = net.dat['node_info'][inst_rc]
    all_cats = [x for x in inst_keys if 'cat-' in x]

    for i in range(len(net.dat['nodes'][inst_rc])):
      # print(net.dat['nodes'][inst_rc])
      inst_dict = {}
      inst_dict['name'] = net.dat['nodes'][inst_rc][i]
      # inst_dict['ini'] = net.dat['node_info'][inst_rc]['ini'][i]
      inst_dict['clust'] = net.dat['node_info'][inst_rc]['clust'].index(i)
      # inst_dict['rank'] = net.dat['node_info'][inst_rc]['rank'][i]
      # if 'rankvar' in inst_keys:
      #   inst_dict['rankvar'] = net.dat['node_info'][inst_rc]['rankvar'][i]
      # fix for similarity matrix
      if len(all_cats) > 0:

        for inst_name_cat in all_cats:

          actual_cat_name = net.dat['node_info'][inst_rc][inst_name_cat][i]
          inst_dict[inst_name_cat] = actual_cat_name

          check_pval = 'pval_'+inst_name_cat.replace('-','_')

          if check_pval in net.dat['node_info'][inst_rc]:
            tmp_pval_name = inst_name_cat.replace('-','_') + '_pval'
            inst_dict[tmp_pval_name] = net.dat['node_info'][inst_rc][check_pval][actual_cat_name]

          tmp_index_name = inst_name_cat.replace('-', '_') + '_index'

          inst_dict[tmp_index_name] = net.dat['node_info'][inst_rc] \
              [tmp_index_name][i]


      if len(net.dat['node_info'][inst_rc]['value']) > 0:
        inst_dict['value'] = net.dat['node_info'][inst_rc]['value'][i]

      if len(net.dat['node_info'][inst_rc]['info']) > 0:
        inst_dict['info'] = net.dat['node_info'][inst_rc]['info'][i]

      if dendro is True:
        inst_dict['group'] = []
        for tmp_dist in all_dist:
          tmp_dist = str(tmp_dist).replace('.', '')
          tmp_append = float(
              net.dat['node_info'][inst_rc]['group'][tmp_dist][i])
          inst_dict['group'].append(tmp_append)

      net.viz[inst_rc + '_nodes'].append(inst_dict)

  # save data as links or mat
  ###########################
  if links is True:
    for i in range(len(net.dat['nodes']['row'])):
      for j in range(len(net.dat['nodes']['col'])):

        inst_dict = {}
        inst_dict['source'] = i
        inst_dict['target'] = j
        inst_dict['value'] = float(net.dat['mat'][i, j])

        if np.isnan(inst_dict['value_orig']):
          inst_dict['value_orig'] = 'NaN'

        net.viz['links'].append(inst_dict)

  else:
    net.viz['mat'] = net.dat['mat'].tolist()




# import numpy as np
# # Assuming calc_clust is in the same directory or accessible
# # from . import calc_clust 

# def viz_json(net, dendro=True, links=False):
#     """
#     Corrected and robustly populates the net.viz dictionary for visualization.

#     This version:
#     1.  Properly populates the pre-initialized net.viz dictionary without overwriting it.
#     2.  Includes checks to diagnose why net.viz['mat'] might be empty.
#     3.  Handles missing data and different data paths gracefully.
#     """
#     # --- 0. Preparation ---
#     # Assumes net.viz has been initialized. We just clear the lists to be repopulated.
#     # This preserves existing keys like 'cat_colors'.

#     # Use hasattr to check for an attribute on an object
#     if not hasattr(net, 'viz') or not isinstance(net.viz, dict):
#         print("🚨 Error: net.viz has not been initialized. Please call initialize_net.viz() first.")
#         net.viz = {} # Initialize as a fallback
    
#     net.viz['row_nodes'] = []
#     net.viz['col_nodes'] = []
#     net.viz['links'] = []
#     net.viz['mat'] = []
#     net.viz['linkage'] = {}

#     # --- 1. Safely Gather Linkage Data ---
#     row_Y = net.dat.get('node_info', {}).get('row', {}).get('Y', np.array([]))
#     col_Y = net.dat.get('node_info', {}).get('col', {}).get('Y', np.array([]))
#     net.viz['linkage'] = {
#         'row': row_Y.tolist(),
#         'col': col_Y.tolist()
#     }

#     # Assuming calc_clust.group_cutoffs() is available
#     # all_dist_cutoffs = calc_clust.group_cutoffs()
#     all_dist_cutoffs = [i / 10.0 for i in range(11)] # Fallback if calc_clust is unavailable

#     # --- 2. Safely Process Nodes (Rows and Columns) ---
#     for inst_rc in ['row', 'col']:
#         nodes = net.dat.get('nodes', {}).get(inst_rc, [])
#         node_info = net.dat.get('node_info', {}).get(inst_rc, {})

#         if not nodes:
#             print(f"⚠️ Warning: No nodes found for '{inst_rc}'. Node list will be empty.")
#             continue

#         # Create safe mappings for sorted orders
#         clust_order = node_info.get('clust', list(range(len(nodes))))
#         rank_order = node_info.get('rank', list(range(len(nodes))))
#         rankvar_order = node_info.get('rankvar', list(range(len(nodes))))
        
#         clust_map = {orig_idx: new_idx for new_idx, orig_idx in enumerate(clust_order)}
#         rank_map = {orig_idx: new_idx for new_idx, orig_idx in enumerate(rank_order)}
#         rankvar_map = {orig_idx: new_idx for new_idx, orig_idx in enumerate(rankvar_order)}

#         for i in range(len(nodes)):
#             inst_dict = {}
#             inst_dict['name'] = nodes[i]
#             # Use .get() with a default list to prevent crashes
#             inst_dict['ini'] = node_info.get('ini', [])[i] if i < len(node_info.get('ini', [])) else -1
#             inst_dict['clust'] = clust_map.get(i, -1)
#             inst_dict['rank'] = rank_map.get(i, -1)
#             if 'rankvar' in node_info:
#                 inst_dict['rankvar'] = rankvar_map.get(i, -1)
            
#             if dendro:
#                 inst_dict['group'] = []
#                 groups = node_info.get('group', {})
#                 for cutoff in all_dist_cutoffs:
#                     key = str(cutoff).replace('.', '')
#                     # Handle keys like '0' -> '00'
#                     if len(key) == 1: key = '0' + key 
#                     group_list = groups.get(key, [])
#                     group_val = group_list[i] if i < len(group_list) else -1.0
#                     inst_dict['group'].append(float(group_val))

#             net.viz[inst_rc + '_nodes'].append(inst_dict)

#     # --- 3. Correctly Populate Matrix or Links ---
#     # Diagnostic Check
#     print(f"DEBUG: 'links' parameter is set to: {links}")
    
#     # CORRECTED: Use the 'in' operator to check for a key in a dictionary-like object.
#     if 'mat' not in net.dat or net.dat['mat'] is None:
#         print("🚨 Error: net.dat['mat'] is missing or None. Cannot generate matrix or links.")
#         return # Exit if there's no data

#     if links is True:
#         print("INFO: Generating data in 'links' format.")
#         source_indices, target_indices = np.where(np.ones_like(net.dat['mat']))
#         values = net.dat['mat'][source_indices, target_indices]
        
#         link_list = [{'source': int(s), 'target': int(t), 'value': v if np.isfinite(v) else None} 
#                      for s, t, v in zip(source_indices, target_indices, values)]
#         net.viz['links'] = link_list
#     # else:
#     #     # This is the default path that populates the matrix.
#     #     print("INFO: Generating data in 'mat' format.")
        
#     #     # Check if the matrix is a numpy array
#     #     if not isinstance(net.dat['mat'], np.ndarray):
#     #         print("⚠️ Warning: net.dat['mat'] is not a numpy array. Attempting to convert.")
#     #         try:
#     #             net.dat['mat'] = np.array(net.dat['mat'])
#     #         except Exception as e:
#     #             print(f"🚨 Error: Could not convert net.dat['mat'] to a numpy array: {e}")
#     #             return
            
#     else:
#         # OPTIMIZED PATH: Send matrix data in a compact, flattened format.
#         print("INFO: Generating data in 'mat' format (OPTIMIZED).")
        
#         # 1. Get the shape of the matrix.
#         mat_shape = net.dat['mat'].shape
        
#         # 2. Flatten the NumPy array into a single, 1D list.
#         #    This is much more memory-efficient than a list of lists.
#         mat_flat = np.where(np.isnan(net.dat['mat']), None, net.dat['mat']).flatten().tolist()

#         # 3. Use the original 'mat' key but with the new, efficient structure.
#         net.viz['mat'] = {
#             'shape': mat_shape, # e.g., [10000, 1000]
#             'flat': mat_flat    # A single long list of values
#         }
        
#         print(f"✅ Successfully generated flattened matrix data for shape: {mat_shape}")


#         # Convert NaN to None for JSON compatibility, then convert to list.
#         mat_for_json = np.where(np.isnan(net.dat['mat']), None, net.dat['mat']).tolist()
#         net.viz['mat'] = mat_for_json
        
#         # Final diagnostic check
#         if not net.viz['mat']:
#             print("🚨 CRITICAL: net.viz['mat'] is still empty after processing. Check the source data in net.dat['mat'].")
#         else:
#             print(f"✅ Successfully generated matrix with dimensions: {len(net.viz['mat'])}x{len(net.viz['mat'][0])}")
