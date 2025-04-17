def cluster_row_and_col(net, dist_type='cosine', linkage_type='average',
                        dendro=True, run_clustering=True, run_rank=True,
                        ignore_cat=False, calc_cat_pval=False, links=False,
                        clust_library='scipy', min_samples=1, min_cluster_size=2):
  ''' cluster net.dat and make visualization json, net.viz.
  optionally leave out dendrogram colorbar groups with dendro argument '''

  # import umap
  import scipy
  from copy import deepcopy
  from scipy.spatial.distance import pdist
  from . import categories, make_viz, cat_pval
  import time
  import traceback
  import numpy as np

  dm = {}
  unique_id = time.time()
  print(f"Function start: cluster_row_and_col called at {unique_id}")
  # traceback.print_stack()

  for axis in ['row', 'col']:

    # save directly to dat structure
    node_info = net.dat['node_info'][axis]

    node_info['ini'] = list(range( len(net.dat['nodes'][axis]), -1, -1))

    tmp_mat = deepcopy(net.dat['mat'])

    # calc distance matrix
    if clust_library != 'hdbscan':

      # dm[axis] = parallel_distance_matrix(tmp_mat, axis, dist_type,n_jobs=4)
      dm[axis] = calc_distance_matrix(tmp_mat, axis, dist_type)

      # Check if NaN exists in dm[axis]
      if np.isnan(dm[axis]).any():
          print("⚠️ Warning: NaN values detected in the distance matrix!")
      else:
          print("✅ No NaN values detected.")


    else:
      dm[axis] = None

    # dm[axis] = calc_distance_matrix(tmp_mat, axis, dist_type)

    # cluster
    if run_clustering is True:
      node_info['clust'], node_info['Y'],  node_info['group'] = clust_and_group(net,
                                                     dm[axis],
                                                     axis,
                                                     tmp_mat,
                                                     dist_type=dist_type,
                                                     linkage_type=linkage_type,
                                                     clust_library=clust_library,
                                                     min_samples=min_samples,
                                                     min_cluster_size=min_cluster_size)
    else:
      dendro = False
      node_info['clust'] = node_info['ini']

    # sorting
    if run_rank is True:
      node_info['rank'] = sort_rank_nodes(net, axis, 'sum')
      node_info['rankvar'] = sort_rank_nodes(net, axis, 'var')
    else:
      node_info['rank'] = node_info['ini']
      node_info['rankvar'] = node_info['ini']

    ##################################
    if ignore_cat is False:
      try:
          categories.calc_cat_clust_order(net, axis)
      except Exception as e:  # Catch all exceptions
          print(f"⚠️ Error in calc_cat_clust_order: {e}")


  if calc_cat_pval is True:
    cat_pval.main(net)

  # make the visualization json
  try:
    make_viz.viz_json(net, dendro, links)
  except Exception as e:  # Catch all exceptions
    print(f"⚠️ Error in calc_cat_clust_order: {e}")


  # make_viz.viz_json(net, dendro, links)

  return dm

# def calc_distance_matrix(tmp_mat, axis, dist_type='cosine'):
#   from scipy.spatial.distance import pdist
#   import numpy as np

#   if axis == 'row':
#     inst_dm = pdist(tmp_mat, metric=dist_type)
#   elif axis == 'col':
#     inst_dm = pdist(tmp_mat.transpose(), metric=dist_type)

#   inst_dm[inst_dm < 0] = float(0)

#   return inst_dm

# def calc_distance_matrix(tmp_mat, axis, dist_type='cosine'):
#     from scipy.spatial.distance import pdist
#     import numpy as np

#     try:
#         if axis not in ['row', 'col']:
#             raise ValueError("Invalid axis. Use 'row' or 'col'.")
        
#         tmp_mat = np.array(tmp_mat, dtype=np.float64)  # Ensure it's a NumPy array

#         print(f"Matrix shape: {tmp_mat.shape}")  

#         if axis == 'row':
#             inst_dm = pdist(tmp_mat, metric=dist_type)
#         elif axis == 'col':
#             inst_dm = pdist(tmp_mat.T, metric=dist_type)  # Using .T instead of transpose() for readability

#         inst_dm[inst_dm < 0] = 0.0  # Ensuring non-negative distances

#         return inst_dm

#     except ValueError as ve:
#         print(f"ValueError: {ve}")
#     except TypeError as te:
#         print(f"TypeError: {te}")
#     except Exception as e:  # Catch-all for unexpected errors
#         print(f"An error occurred: {e}")

#     return None  # Return None in case of failure


def calc_distance_matrix(tmp_mat, axis, dist_type='cosine'):
    from scipy.spatial.distance import pdist
    import numpy as np

    try:
        if axis not in ['row', 'col']:
            raise ValueError("Invalid axis. Use 'row' or 'col'.")

        tmp_mat = np.array(tmp_mat, dtype=np.float64)  # Convert to float (ensures numerical stability)



        # Normalize for cosine distance (prevents division by zero issues)
        if dist_type == 'cosine':
            norms = np.linalg.norm(tmp_mat, axis=1 if axis == 'row' else 0, keepdims=True)
            norms[norms == 0] = 1  # Avoid division by zero
            tmp_mat = tmp_mat / norms  # Normalize

        # Compute distance matrix
        inst_dm = pdist(tmp_mat if axis == 'row' else tmp_mat.T, metric=dist_type)

        # Ensure non-negative distances
        inst_dm[inst_dm < 0] = 0.0  

        return inst_dm

    except ValueError as ve:
        print(f"ValueError: {ve}")
    except TypeError as te:
        print(f"TypeError: {te}")
    except Exception as e:  # Catch-all for unexpected errors
        print(f"An error occurred: {e}")

    return None  # Return None in case of failure




def calculate_blockwise_distances(block_i, block_j, metric='cosine'):
    from scipy.spatial.distance import pdist, cdist, squareform

    """Calculate distances between rows of block_i and block_j."""
    if block_i is block_j:
        # Diagonal block: Calculate the distance within the block
        return squareform(pdist(block_i, metric=metric))
    else:
        # Off-diagonal block: Calculate the distance between block_i and block_j using cdist
        return cdist(block_i, block_j, metric=metric)

def parallel_distance_matrix(data, axis='row', metric='cosine', n_jobs=-1):
    import numpy as np
    from scipy.spatial.distance import squareform
    from joblib import Parallel, delayed
    """Calculate the full pairwise distance matrix in parallel using block-wise processing."""
    
    if axis == 'col':
        # Transpose the matrix to work with columns as rows
        data = data.T
    
    n = len(data)
    # n_jobs = calculate_free_cores()
    print('**** n jobs are as follows ****',n_jobs)

    block_size = n // n_jobs if n_jobs > 1 else n
    blocks = [data[i:i + block_size] for i in range(0, n, block_size)]
    
    # Calculate the distances for each block combination in parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(calculate_blockwise_distances)(blocks[i], blocks[j], metric=metric)
        for i in range(len(blocks)) for j in range(i, len(blocks))
    )
    
    # Combine the results into a full distance matrix
    distance_matrix = np.zeros((n, n))
    
    idx = 0
    for i in range(len(blocks)):
        for j in range(i, len(blocks)):
            block = results[idx]
            if i == j:
                distance_matrix[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size] = block
            else:
                # Handle off-diagonal blocks
                distance_matrix[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size] = block
                distance_matrix[j*block_size:(j+1)*block_size, i*block_size:(i+1)*block_size] = block.T
            idx += 1

    # Convert the full distance matrix to condensed form
    condensed_distance_matrix = squareform(distance_matrix, checks=False)

    return condensed_distance_matrix

def clust_and_group(net, inst_dm, axis, mat, dist_type='cosine', linkage_type='average',
                    clust_library='scipy', min_samples=1, min_cluster_size=2):

  # print(clust_library)

  import scipy.cluster.hierarchy as hier
  import pandas as pd

  ### Added extra for debugging delete it #####
  import numpy as np
  if clust_library == 'scipy':
    # Check for infinite values
    # has_nan = np.isnan(inst_dm).any()
    # if has_nan:
    #     print("NaN values found in the matrix.")
    #     nan_indices = np.where(np.isnan(inst_dm))
    #     print("Indices of NaN values:", nan_indices)
    # else:
    #     print("No NaN values found in the matrix.")


    Y = hier.linkage(inst_dm, method=linkage_type)

  elif clust_library == 'fastcluster':
    import fastcluster
    Y = fastcluster.linkage(inst_dm, method=linkage_type)

  elif clust_library == 'hdbscan':
    # print('HDBSCAN!')
    import hdbscan


    # pca-umap-hdbscan using data (no pre-cal distance matrix)
    ######################################################
    from sklearn.decomposition import PCA
    clusterer = hdbscan.HDBSCAN(min_samples=min_samples,
                                min_cluster_size=min_cluster_size)

    # rows are the data points, cols are dimensions
    n_components = 50
    if axis == 'row':
      if mat.shape[1] > n_components:
        low_d_mat = PCA(n_components=n_components).fit_transform(mat)
      else:
        low_d_mat = mat

    elif axis == 'col':
      if mat.shape[0] > n_components:
        low_d_mat = PCA(n_components=n_components).fit_transform(mat.transpose())
      else:
        low_d_mat = mat.transpose()


    # run UMAP on low_d_mat (after PCA)
    # print('running umap!!!!!!!!!!!!!!!!!!!!!!!!!!')
    import umap
    umap_mat = umap.UMAP(
                          metric=dist_type,
                          n_neighbors=5,
                          min_dist=0.0,
                          n_components=2,
                          random_state=42,
                          ).fit_transform(low_d_mat)

    umap_df = pd.DataFrame(umap_mat.transpose(),
                           index=['x','y'],
                           columns=net.dat['nodes'][axis])

    net.umap[axis] = umap_df
    clusterer.fit(umap_mat)

    Y = clusterer.single_linkage_tree_.to_numpy()
    Z = hier.dendrogram(Y, no_plot=True)

  Z = hier.dendrogram(Y, no_plot=True)
  # Z = hier.dendrogram(Y, no_plot=True, truncate_mode='level', p=12)

  # if axis == 'row':
  #   print(Z)

  inst_clust_order = Z['leaves']
  # Only calculate groups for hierarchical clustering methods
  groups = {}
  if clust_library in ['scipy', 'fastcluster']:
      all_dist = group_cutoffs()
      for inst_dist in all_dist:
          inst_key = str(inst_dist).replace('.', '')
          groups[inst_key] = hier.fcluster(Y, inst_dist * inst_dm.max(), 'distance')
          groups[inst_key] = groups[inst_key].tolist()
  elif clust_library == 'hdbscan':
      # For HDBSCAN, you might want to handle groups differently
      # One option is to use the cluster labels at different probabilities
      # This is just a placeholder - you'll need to adapt this for HDBSCAN
      groups['hdbscan_clusters'] = clusterer.labels_.tolist()

  return inst_clust_order, Y, groups

def sort_rank_nodes(net, rowcol, rank_type):
  import numpy as np
  from operator import itemgetter
  from copy import deepcopy

  tmp_nodes = deepcopy(net.dat['nodes'][rowcol])
  inst_mat = deepcopy(net.dat['mat'])

  sum_term = []
  for i in range(len(tmp_nodes)):
    inst_dict = {}
    inst_dict['name'] = tmp_nodes[i]

    if rowcol == 'row':
      if rank_type == 'sum':
        inst_dict['rank'] = np.sum(inst_mat[i, :])
      elif rank_type == 'var':
        inst_dict['rank'] = np.var(inst_mat[i, :])
    else:
      if rank_type == 'sum':
        inst_dict['rank'] = np.sum(inst_mat[:, i])
      elif rank_type == 'var':
        inst_dict['rank'] = np.var(inst_mat[:, i])

    sum_term.append(inst_dict)

  sum_term = sorted(sum_term, key=itemgetter('rank'), reverse=False)

  tmp_sort_nodes = []
  for inst_dict in sum_term:
    tmp_sort_nodes.append(inst_dict['name'])

  sort_index = []
  for inst_node in tmp_nodes:
    sort_index.append(tmp_sort_nodes.index(inst_node))

  return sort_index

def group_cutoffs():
  all_dist = []
  for i in range(11):
    all_dist.append(float(i) / 10)
  return all_dist
