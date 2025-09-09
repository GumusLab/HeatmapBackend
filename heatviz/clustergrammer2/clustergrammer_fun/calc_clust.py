# def cluster_row_and_col(net, dist_type='cosine', linkage_type='average',
#                         dendro=True, run_clustering=True, run_rank=True,
#                         ignore_cat=False, calc_cat_pval=False, links=False,
#                         clust_library='scipy', min_samples=1, min_cluster_size=2):
#   ''' cluster net.dat and make visualization json, net.viz.
#   optionally leave out dendrogram colorbar groups with dendro argument '''

#   # import umap
#   import scipy
#   from copy import deepcopy
#   from scipy.spatial.distance import pdist
#   from . import categories, make_viz, cat_pval
#   import time
#   import traceback
#   import numpy as np

#   dm = {}
#   unique_id = time.time()
#   print(f"Function start: cluster_row_and_col called at {unique_id}")
#   # traceback.print_stack()

#   for axis in ['row', 'col']:

#     # save directly to dat structure
#     node_info = net.dat['node_info'][axis]

#     node_info['ini'] = list(range( len(net.dat['nodes'][axis]), -1, -1))

#     tmp_mat = deepcopy(net.dat['mat'])

#     # calc distance matrix
#     if clust_library != 'hdbscan':

#       dm[axis] = parallel_distance_matrix(tmp_mat, axis, dist_type,n_jobs=4)
#       # dm[axis] = calc_distance_matrix(tmp_mat, axis, dist_type)

#       # Check if NaN exists in dm[axis]
#       if np.isnan(dm[axis]).any():
#           print("⚠️ Warning: NaN values detected in the distance matrix!")
#       else:
#           print("✅ No NaN values detected.")


#     else:
#       dm[axis] = None

#     # dm[axis] = calc_distance_matrix(tmp_mat, axis, dist_type)

#     # cluster
#     if run_clustering is True:
#       node_info['clust'], node_info['Y'],  node_info['group'] = clust_and_group(net,
#                                                      dm[axis],
#                                                      axis,
#                                                      tmp_mat,
#                                                      dist_type=dist_type,
#                                                      linkage_type=linkage_type,
#                                                      clust_library=clust_library,
#                                                      min_samples=min_samples,
#                                                      min_cluster_size=min_cluster_size)
#     else:
#       dendro = False
#       node_info['clust'] = node_info['ini']

#     # sorting
#     if run_rank is True:
#       node_info['rank'] = sort_rank_nodes(net, axis, 'sum')
#       node_info['rankvar'] = sort_rank_nodes(net, axis, 'var')
#     else:
#       node_info['rank'] = node_info['ini']
#       node_info['rankvar'] = node_info['ini']

#     ##################################
#     if ignore_cat is False:
#       try:
#           categories.calc_cat_clust_order(net, axis)
#       except Exception as e:  # Catch all exceptions
#           print(f"⚠️ Error in calc_cat_clust_order: {e}")


#   if calc_cat_pval is True:
#     cat_pval.main(net)

#   # make the visualization json
#   try:
#     make_viz.viz_json(net, dendro, links)
#   except Exception as e:  # Catch all exceptions
#     print(f"⚠️ Error in calc_cat_clust_order: {e}")


#   # make_viz.viz_json(net, dendro, links)

#   return dm



# # def calc_distance_matrix(tmp_mat, axis, dist_type='cosine'):
# #     from scipy.spatial.distance import pdist
# #     import numpy as np

# #     try:
# #         if axis not in ['row', 'col']:
# #             raise ValueError("Invalid axis. Use 'row' or 'col'.")

# #         tmp_mat = np.array(tmp_mat, dtype=np.float64)  # Convert to float (ensures numerical stability)



# #         # Normalize for cosine distance (prevents division by zero issues)
# #         if dist_type == 'cosine':
# #             norms = np.linalg.norm(tmp_mat, axis=1 if axis == 'row' else 0, keepdims=True)
# #             norms[norms == 0] = 1  # Avoid division by zero
# #             tmp_mat = tmp_mat / norms  # Normalize

# #         # Compute distance matrix
# #         inst_dm = pdist(tmp_mat if axis == 'row' else tmp_mat.T, metric=dist_type)

# #         # Ensure non-negative distances
# #         inst_dm[inst_dm < 0] = 0.0  

# #         return inst_dm

# #     except ValueError as ve:
# #         print(f"ValueError: {ve}")
# #     except TypeError as te:
# #         print(f"TypeError: {te}")
# #     except Exception as e:  # Catch-all for unexpected errors
# #         print(f"An error occurred: {e}")

# #     return None  # Return None in case of failure

# def calc_distance_matrix(tmp_mat, axis, dist_type='cosine'):
#     from scipy.spatial.distance import pdist
#     import numpy as np

#     try:
#         if axis not in ['row', 'col']:
#             raise ValueError("Invalid axis. Use 'row' or 'col'.")

#         # Convert to numpy array
#         tmp_mat = np.array(tmp_mat, dtype=np.float64)
        
#         # Check for and handle NaN values in the input matrix
#         if np.isnan(tmp_mat).any():
#             print("⚠️ NaN values found in input matrix. Replacing with zeros...")
#             tmp_mat = np.nan_to_num(tmp_mat, nan=0.0)
        
#         # Check for zero rows/columns which cause problems with cosine distance
#         if axis == 'row':
#             # Calculate row norms
#             row_norms = np.linalg.norm(tmp_mat, axis=1)
#             zero_rows = np.where(row_norms == 0)[0]
#             if len(zero_rows) > 0:
#                 print(f"⚠️ Found {len(zero_rows)} zero-norm rows which will cause NaN in cosine distance")
#                 # Add small random noise to zero rows to prevent NaN
#                 for idx in zero_rows:
#                     tmp_mat[idx] = np.random.normal(0, 1e-6, tmp_mat.shape[1])
                
#             # Compute distance matrix
#             inst_dm = pdist(tmp_mat, metric=dist_type)
#         elif axis == 'col':
#             # Calculate column norms
#             col_norms = np.linalg.norm(tmp_mat, axis=0)
#             zero_cols = np.where(col_norms == 0)[0]
#             if len(zero_cols) > 0:
#                 print(f"⚠️ Found {len(zero_cols)} zero-norm columns which will cause NaN in cosine distance")
#                 # Add small random noise to zero columns to prevent NaN
#                 for idx in zero_cols:
#                     tmp_mat[:, idx] = np.random.normal(0, 1e-6, tmp_mat.shape[0])
                
#             # Compute distance matrix
#             inst_dm = pdist(tmp_mat.T, metric=dist_type)
        
#         # Check for NaN or Inf values in the distance matrix and fix them
#         if not np.isfinite(inst_dm).all():
#             print(f"⚠️ Found {np.sum(~np.isfinite(inst_dm))} non-finite values in distance matrix. Replacing them...")
#             inst_dm = np.nan_to_num(inst_dm, nan=0.0, posinf=1.0, neginf=0.0)
        
#         # Ensure non-negative distances
#         inst_dm[inst_dm < 0] = 0.0

#         return inst_dm

#     except Exception as e:
#         print(f"❌ Error in calc_distance_matrix: {e}")
#         import traceback
#         traceback.print_exc()
        
#         # Return a fake distance matrix as fallback (all zeros except diagonal)
#         size = tmp_mat.shape[0] if axis == 'row' else tmp_mat.shape[1]
#         n_distances = size * (size - 1) // 2  # Size of condensed distance matrix
#         return np.zeros(n_distances)




# # def calculate_blockwise_distances(block_i, block_j, metric='cosine'):
# #     from scipy.spatial.distance import pdist, cdist, squareform

# #     """Calculate distances between rows of block_i and block_j."""
# #     if block_i is block_j:
# #         # Diagonal block: Calculate the distance within the block
# #         return squareform(pdist(block_i, metric=metric))
# #     else:
# #         # Off-diagonal block: Calculate the distance between block_i and block_j using cdist
# #         return cdist(block_i, block_j, metric=metric)
    

# def calculate_blockwise_distances(block_i, block_j, metric='cosine'):
#     from scipy.spatial.distance import pdist, cdist, squareform
#     import numpy as np
    
#     try:
#         # Check inputs for NaN values
#         block_i = np.nan_to_num(block_i, nan=0.0)
#         block_j = np.nan_to_num(block_j, nan=0.0)
        
#         # Special handling for cosine with zero vectors
#         if metric == 'cosine':
#             # Check and fix zero norm vectors in block_i
#             norms_i = np.linalg.norm(block_i, axis=1)
#             zero_i = np.where(norms_i < 1e-10)[0]
#             if len(zero_i) > 0:
#                 for idx in zero_i:
#                     block_i[idx] = np.random.normal(0, 1e-6, block_i.shape[1])
            
#             # If comparing different blocks, check block_j too
#             if block_i is not block_j:
#                 norms_j = np.linalg.norm(block_j, axis=1)
#                 zero_j = np.where(norms_j < 1e-10)[0]
#                 if len(zero_j) > 0:
#                     for idx in zero_j:
#                         block_j[idx] = np.random.normal(0, 1e-6, block_j.shape[1])
        
#         # Calculate distances
#         if block_i is block_j:
#             # For same block, use pdist (more efficient)
#             result = squareform(pdist(block_i, metric=metric))
#         else:
#             # For different blocks, use cdist
#             result = cdist(block_i, block_j, metric=metric)
        
#         # Check result for NaN/Inf values
#         if not np.isfinite(result).all():
#             print(f"⚠️ Non-finite values in block result. Fixing...")
#             result = np.nan_to_num(result, nan=0.0, posinf=1.0, neginf=0.0)
        
#         return result
        
#     except Exception as e:
#         print(f"❌ Error in calculate_blockwise_distances: {e}")
#         # Return zeros as fallback
#         shape = (len(block_i), len(block_j)) if block_i is not block_j else (len(block_i), len(block_i))
#         return np.zeros(shape)

# # def parallel_distance_matrix(data, axis='row', metric='cosine', n_jobs=-1):
# #     import numpy as np
# #     from scipy.spatial.distance import squareform
# #     from joblib import Parallel, delayed
# #     """Calculate the full pairwise distance matrix in parallel using block-wise processing."""
    
# #     if axis == 'col':
# #         # Transpose the matrix to work with columns as rows
# #         data = data.T
    
# #     n = len(data)
# #     # n_jobs = calculate_free_cores()
# #     print('**** n jobs are as follows ****',n_jobs)

# #     block_size = n // n_jobs if n_jobs > 1 else n
# #     blocks = [data[i:i + block_size] for i in range(0, n, block_size)]
    
# #     # Calculate the distances for each block combination in parallel
# #     results = Parallel(n_jobs=n_jobs)(
# #         delayed(calculate_blockwise_distances)(blocks[i], blocks[j], metric=metric)
# #         for i in range(len(blocks)) for j in range(i, len(blocks))
# #     )
    
# #     # Combine the results into a full distance matrix
# #     distance_matrix = np.zeros((n, n))
    
# #     idx = 0
# #     for i in range(len(blocks)):
# #         for j in range(i, len(blocks)):
# #             block = results[idx]
# #             if i == j:
# #                 distance_matrix[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size] = block
# #             else:
# #                 # Handle off-diagonal blocks
# #                 distance_matrix[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size] = block
# #                 distance_matrix[j*block_size:(j+1)*block_size, i*block_size:(i+1)*block_size] = block.T
# #             idx += 1

# #     # Convert the full distance matrix to condensed form
# #     condensed_distance_matrix = squareform(distance_matrix, checks=False)

# #     return condensed_distance_matrix


# def parallel_distance_matrix(data, axis='row', metric='cosine', n_jobs=4):
#     import numpy as np
#     from scipy.spatial.distance import squareform, pdist
#     from joblib import Parallel, delayed
    
#     try:
#         # Handle the case if data is not a numpy array
#         data = np.array(data, dtype=np.float64)
        
#         # Check for and fix NaN values
#         if np.isnan(data).any():
#             print(f"⚠️ Found {np.sum(np.isnan(data))} NaN values in input data. Replacing with zeros...")
#             data = np.nan_to_num(data, nan=0.0)
        
#         # Transpose if working with columns
#         if axis == 'col':
#             data = data.T
        
#         # Handle zero-norm vectors which cause issues with cosine distance
#         if metric == 'cosine':
#             norms = np.linalg.norm(data, axis=1)
#             zero_indices = np.where(norms < 1e-10)[0]
#             if len(zero_indices) > 0:
#                 print(f"⚠️ Found {len(zero_indices)} zero-norm vectors. Adding small noise...")
#                 for idx in zero_indices:
#                     data[idx] = np.random.normal(0, 1e-6, data.shape[1])
        
#         # For smaller matrices, use regular pdist for simplicity and reliability
#         n = len(data)
#         if n < 1000 or n_jobs <= 1:
#             condensed = pdist(data, metric=metric)
#             # Check for and fix non-finite values
#             if not np.isfinite(condensed).all():
#                 print(f"⚠️ Found {np.sum(~np.isfinite(condensed))} non-finite values. Fixing...")
#                 condensed = np.nan_to_num(condensed, nan=0.0, posinf=1.0, neginf=0.0)
#             return condensed
        
#         # For larger matrices, use block-wise parallel computation
#         block_size = max(1, n // n_jobs)
#         blocks = [data[i:min(i + block_size, n)] for i in range(0, n, block_size)]
        
#         # Calculate blocks in parallel
#         results = Parallel(n_jobs=n_jobs)(
#             delayed(calculate_blockwise_distances)(blocks[i], blocks[j], metric=metric)
#             for i in range(len(blocks)) for j in range(i, len(blocks))
#         )
        
#         # Build the full matrix
#         full_matrix = np.zeros((n, n))
#         idx = 0
#         for i in range(len(blocks)):
#             block_i_start = i * block_size
#             block_i_end = min(block_i_start + len(blocks[i]), n)
            
#             for j in range(i, len(blocks)):
#                 block_j_start = j * block_size
#                 block_j_end = min(block_j_start + len(blocks[j]), n)
                
#                 block_result = results[idx]
                
#                 # Check for NaN values in this block
#                 if np.isnan(block_result).any():
#                     print(f"⚠️ NaN values in block {i},{j}. Fixing...")
#                     block_result = np.nan_to_num(block_result, nan=0.0)
                
#                 # Set the block in the full matrix
#                 full_matrix[block_i_start:block_i_end, block_j_start:block_j_end] = block_result
                
#                 # Set the symmetric block
#                 if i != j:
#                     full_matrix[block_j_start:block_j_end, block_i_start:block_i_end] = block_result.T
                
#                 idx += 1
        
#         # Convert to condensed form
#         condensed = squareform(full_matrix, checks=False)
        
#         # Final check for bad values
#         if not np.isfinite(condensed).all():
#             print(f"⚠️ Found {np.sum(~np.isfinite(condensed))} non-finite values in final result. Fixing...")
#             condensed = np.nan_to_num(condensed, nan=0.0, posinf=1.0, neginf=0.0)
        
#         return condensed
        
#     except Exception as e:
#         print(f"❌ Error in parallel_distance_matrix: {e}")
#         import traceback
#         traceback.print_exc()
        
#         # Return a fallback distance matrix
#         n = len(data) if isinstance(data, np.ndarray) else 0
#         if n > 0:
#             n_distances = n * (n - 1) // 2
#             return np.zeros(n_distances)
#         else:
#             return np.array([])

# def clust_and_group(net, inst_dm, axis, mat, dist_type='cosine', linkage_type='average',
#                     clust_library='scipy', min_samples=1, min_cluster_size=2):

#   # print(clust_library)

#   import scipy.cluster.hierarchy as hier
#   import pandas as pd

#   ### Added extra for debugging delete it #####
#   import numpy as np
#   if clust_library == 'scipy':
#     # Check for infinite values
#     # has_nan = np.isnan(inst_dm).any()
#     # if has_nan:
#     #     print("NaN values found in the matrix.")
#     #     nan_indices = np.where(np.isnan(inst_dm))
#     #     print("Indices of NaN values:", nan_indices)
#     # else:
#     #     print("No NaN values found in the matrix.")


#     Y = hier.linkage(inst_dm, method=linkage_type)

#   elif clust_library == 'fastcluster':
#     import fastcluster
#     Y = fastcluster.linkage(inst_dm, method=linkage_type)

#   elif clust_library == 'hdbscan':
#     # print('HDBSCAN!')
#     import hdbscan


#     # pca-umap-hdbscan using data (no pre-cal distance matrix)
#     ######################################################
#     from sklearn.decomposition import PCA
#     clusterer = hdbscan.HDBSCAN(min_samples=min_samples,
#                                 min_cluster_size=min_cluster_size)

#     # rows are the data points, cols are dimensions
#     n_components = 50
#     if axis == 'row':
#       if mat.shape[1] > n_components:
#         low_d_mat = PCA(n_components=n_components).fit_transform(mat)
#       else:
#         low_d_mat = mat

#     elif axis == 'col':
#       if mat.shape[0] > n_components:
#         low_d_mat = PCA(n_components=n_components).fit_transform(mat.transpose())
#       else:
#         low_d_mat = mat.transpose()


#     # run UMAP on low_d_mat (after PCA)
#     # print('running umap!!!!!!!!!!!!!!!!!!!!!!!!!!')
#     import umap
#     umap_mat = umap.UMAP(
#                           metric=dist_type,
#                           n_neighbors=5,
#                           min_dist=0.0,
#                           n_components=2,
#                           random_state=42,
#                           ).fit_transform(low_d_mat)

#     umap_df = pd.DataFrame(umap_mat.transpose(),
#                            index=['x','y'],
#                            columns=net.dat['nodes'][axis])

#     net.umap[axis] = umap_df
#     clusterer.fit(umap_mat)

#     Y = clusterer.single_linkage_tree_.to_numpy()
#     Z = hier.dendrogram(Y, no_plot=True)

#   Z = hier.dendrogram(Y, no_plot=True)
#   # Z = hier.dendrogram(Y, no_plot=True, truncate_mode='level', p=12)

#   # if axis == 'row':
#   #   print(Z)

#   inst_clust_order = Z['leaves']
#   # Only calculate groups for hierarchical clustering methods
#   groups = {}
#   if clust_library in ['scipy', 'fastcluster']:
#       all_dist = group_cutoffs()
#       for inst_dist in all_dist:
#           inst_key = str(inst_dist).replace('.', '')
#           groups[inst_key] = hier.fcluster(Y, inst_dist * inst_dm.max(), 'distance')
#           groups[inst_key] = groups[inst_key].tolist()
#   elif clust_library == 'hdbscan':
#       # For HDBSCAN, you might want to handle groups differently
#       # One option is to use the cluster labels at different probabilities
#       # This is just a placeholder - you'll need to adapt this for HDBSCAN
#       groups['hdbscan_clusters'] = clusterer.labels_.tolist()

#   return inst_clust_order, Y, groups

# def sort_rank_nodes(net, rowcol, rank_type):
#   import numpy as np
#   from operator import itemgetter
#   from copy import deepcopy

#   tmp_nodes = deepcopy(net.dat['nodes'][rowcol])
#   inst_mat = deepcopy(net.dat['mat'])

#   sum_term = []
#   for i in range(len(tmp_nodes)):
#     inst_dict = {}
#     inst_dict['name'] = tmp_nodes[i]

#     if rowcol == 'row':
#       if rank_type == 'sum':
#         inst_dict['rank'] = np.sum(inst_mat[i, :])
#       elif rank_type == 'var':
#         inst_dict['rank'] = np.var(inst_mat[i, :])
#     else:
#       if rank_type == 'sum':
#         inst_dict['rank'] = np.sum(inst_mat[:, i])
#       elif rank_type == 'var':
#         inst_dict['rank'] = np.var(inst_mat[:, i])

#     sum_term.append(inst_dict)

#   sum_term = sorted(sum_term, key=itemgetter('rank'), reverse=False)

#   tmp_sort_nodes = []
#   for inst_dict in sum_term:
#     tmp_sort_nodes.append(inst_dict['name'])

#   sort_index = []
#   for inst_node in tmp_nodes:
#     sort_index.append(tmp_sort_nodes.index(inst_node))

#   return sort_index

# def group_cutoffs():
#   all_dist = []
#   for i in range(11):
#     all_dist.append(float(i) / 10)
#   return all_dist


import numpy as np
import scipy.cluster.hierarchy as hier
from scipy.spatial.distance import pdist, squareform
from joblib import Parallel, delayed, effective_n_jobs
from copy import deepcopy
import time
import traceback
import pandas as pd
from .gene_imputation import GeneExpressionImputer  # Add this line



# def cluster_row_and_col(net, dist_type='cosine', linkage_type='average',
#                         dendro=True, run_clustering=True, run_rank=True,
#                         ignore_cat=False, calc_cat_pval=False, links=False,
#                         clust_library='scipy', min_samples=1, min_cluster_size=2,
#                         n_jobs=-1, max_cache_size_mb=500):
def cluster_row_and_col(net, dist_type='cosine', linkage_type='average',
                        dendro=True, run_clustering=True, run_rank=True,
                        ignore_cat=False, calc_cat_pval=False, links=False,
                        clust_library='scipy', min_samples=1, min_cluster_size=2,
                        n_jobs=-1, max_cache_size_mb=500,
                        imputation_method='auto', **imputation_kwargs):
    """
    Optimized clustering function with better memory management and performance.
    
    Parameters:
    -----------
    max_cache_size_mb : int, default 500
        Maximum memory to use for caching distance matrices (in MB). 
        Set to 0 to disable caching. For large matrices (>10K rows), 
        caching may be disabled automatically.
    n_jobs : int, default -1
        Number of parallel jobs. -1 uses all available cores
    """
    
    from . import categories, make_viz, cat_pval
    
    dm = {}
    distance_cache = {}
    unique_id = time.time()
    print(f"Function start: cluster_row_and_col called at {unique_id}")
    
    # Pre-process matrix once to avoid repeated deep copies
    base_mat = np.array(net.dat['mat'], dtype=np.float64)
    # Check for any NaN values before imputation
    nan_count = np.isnan(base_mat).sum()
    if nan_count > 0:
        nan_rows = np.isnan(base_mat).any(axis=1).sum()
        nan_cols = np.isnan(base_mat).any(axis=0).sum()
        print(f'Rows with NaN: {nan_rows}')
        print(f'Cols with NaN: {nan_cols}')

    # base_mat = _preprocess_matrix(base_mat)
    base_mat = _preprocess_matrix(base_mat, imputation_method, **imputation_kwargs)

    
    # Check if shapes match the node info
    expected_rows = len(net.dat["node_info"]["row"]["full_names"])
    expected_cols = len(net.dat["node_info"]["col"]["full_names"])
    
    
    if base_mat.shape[0] != expected_rows:
        print(f'ROW MISMATCH: Matrix has {base_mat.shape[0]} rows but node_info has {expected_rows}')
    else:
        print(f'Row count matches')
        
    if base_mat.shape[1] != expected_cols:
        print(f'COL MISMATCH: Matrix has {base_mat.shape[1]} cols but node_info has {expected_cols}')
    else:
        print(f'Column count matches')
    
    # Smart caching decision based on matrix size
    use_caching = _should_cache_distances(base_mat, max_cache_size_mb)
    if not use_caching:
        print(f"⚠️ Matrix too large ({base_mat.shape}) - disabling distance matrix caching to save memory")
    
    for axis in ['row', 'col']:
        print(f"Processing {axis} clustering...")
        
        # Initialize node info
        node_info = net.dat['node_info'][axis]
        node_info['ini'] = list(range(len(net.dat['nodes'][axis]) - 1, -1, -1))
        
        # Calculate or retrieve cached distance matrix
        if clust_library != 'hdbscan':
            cache_key = f"{axis}_{dist_type}" if use_caching else None
            
            if use_caching and cache_key in distance_cache:
                print(f"Using cached distance matrix for {axis}")
                dm[axis] = distance_cache[cache_key]
            else:
                dm[axis] = optimized_distance_matrix(
                    base_mat, axis, dist_type, n_jobs=n_jobs
                )
                if use_caching and cache_key:
                    distance_cache[cache_key] = dm[axis]
            
            # Validate distance matrix
            if not _validate_distance_matrix(dm[axis]):
                print(f"⚠️ Invalid distance matrix for {axis}, using fallback")
                dm[axis] = _create_fallback_distances(base_mat, axis)
        else:
            dm[axis] = None
        
        # Clustering
        if run_clustering:
            node_info['clust'], node_info['Y'], node_info['group'] = clust_and_group(
                net, dm[axis], axis, base_mat,
                dist_type=dist_type, linkage_type=linkage_type,
                clust_library=clust_library, min_samples=min_samples,
                min_cluster_size=min_cluster_size
            )
        else:
            dendro = False
            node_info['clust'] = node_info['ini']
        
        # Ranking - optimized to avoid repeated matrix operations
        if run_rank:
            node_info['rank'] = optimized_sort_rank_nodes(base_mat, axis, 'sum')
            node_info['rankvar'] = optimized_sort_rank_nodes(base_mat, axis, 'var')
        else:
            node_info['rank'] = node_info['ini']
            node_info['rankvar'] = node_info['ini']
        
        # Category calculations
        if not ignore_cat:
            try:
                # Check if this axis has meaningful categories before processing
                inst_keys = net.dat['node_info'][axis]
                all_cats = [x for x in inst_keys if 'cat-' in x]
                
                has_meaningful_cats = False
                for cat_name in all_cats:
                    dict_name = 'dict_' + cat_name.replace('-', '_')
                    if dict_name in inst_keys:
                        dict_cat = inst_keys[dict_name]
                        if isinstance(dict_cat, dict):
                            # Check if all categories are NaN
                            nan_count = sum(1 for cat in dict_cat.keys() if str(cat) == 'nan')
                            if nan_count < len(dict_cat):
                                has_meaningful_cats = True
                                break
                
                if has_meaningful_cats:
                    categories.calc_cat_clust_order(net, axis)
                else:
                    print(f"🔍 DEBUG: Skipping category processing for {axis} - all categories are NaN")
                    
            except Exception as e:
                print(f"⚠️ Error in calc_cat_clust_order: {e}")
    
    # Category p-values
    if calc_cat_pval:
        try:
            cat_pval.main(net)
        except Exception as e:
            print(f"⚠️ Error in cat_pval: {e}")
    
    # Visualization
    try:
        make_viz.viz_json(net, dendro, links)
    except Exception as e:
        print(f"⚠️ Error in make_viz: {e}")

    print(dm)
    
    return dm


def _should_cache_distances(mat, max_cache_mb):
    """
    Determine if distance matrix caching is advisable based on matrix size and memory constraints.
    
    For a matrix with n rows/cols, distance matrix uses n*(n-1)/2 * 8 bytes (float64).
    """
    n_rows, n_cols = mat.shape
    
    # Calculate memory requirements for both row and column distance matrices
    row_distances = n_rows * (n_rows - 1) // 2
    col_distances = n_cols * (n_cols - 1) // 2
    
    # Memory in MB (8 bytes per float64)
    total_cache_mb = (row_distances + col_distances) * 8 / (1024 * 1024)
    
    print(f"Estimated cache memory: {total_cache_mb:.1f} MB for matrix {mat.shape}")
    
    # Don't cache if it would exceed the limit
    if total_cache_mb > max_cache_mb:
        return False
    
    # Also don't cache for very large matrices (>15K in any dimension)
    if max(n_rows, n_cols) > 15000:
        return False
    
    return True


# def _preprocess_matrix(mat):
#     """Preprocess matrix to handle common issues that cause clustering problems."""
#     # Handle NaN values
#     if np.isnan(mat).any():
#         print(f"⚠️ Found {np.sum(np.isnan(mat))} NaN values. Replacing with zeros...")
#         mat = np.nan_to_num(mat, nan=0.0, posinf=1.0, neginf=0.0)
    
#     # Handle infinite values
#     if not np.isfinite(mat).all():
#         print("⚠️ Found infinite values. Clipping to finite range...")
#         mat = np.clip(mat, -1e10, 1e10)
    
#     return mat

def _preprocess_matrix(mat, imputation_method='auto', **imputation_kwargs):
    """Preprocess matrix with sophisticated imputation for gene expression data."""
    # Handle infinite values first
    if not np.isfinite(mat).all():
        print("⚠️ Found infinite values. Clipping to finite range...")
        mat[~np.isfinite(mat)] = np.nan

        # mat = np.nan_to_num(mat, nan=np.nan, posinf=1e10, neginf=-1e10)
    
    # Check if we have missing values
    if np.isnan(mat).any():
        print(f"⚠️ Found {np.sum(np.isnan(mat))} missing values. Performing {imputation_method} imputation...")
        
        try:
            # Initialize imputer
            imputer = GeneExpressionImputer(method=imputation_method, **imputation_kwargs)
            
            # Perform imputation
            mat = imputer.fit_transform(mat)
            
            # Verify no missing values remain
            if np.isnan(mat).any():
                print("⚠️ Imputation failed to remove all missing values. Using fallback...")
                mat = np.nan_to_num(mat, nan=0.0)
                
        except Exception as e:
            print(f"⚠️ Imputation failed with error: {e}. Using simple replacement...")
            mat = np.nan_to_num(mat, nan=0.0, posinf=1.0, neginf=0.0)

    print('********** Imputation is done *************')
    
    return mat


def _validate_distance_matrix(dm):
    """Validate that distance matrix is reasonable."""
    if dm is None or len(dm) == 0:
        return False
    
    return (np.isfinite(dm).all() and 
            (dm >= 0).all() and 
            not np.isnan(dm).any())


def _create_fallback_distances(mat, axis):
    """Create a simple fallback distance matrix."""
    n = mat.shape[0] if axis == 'row' else mat.shape[1]
    n_distances = n * (n - 1) // 2
    return np.random.uniform(0, 1, n_distances)  # Random but valid distances


def optimized_distance_matrix(data, axis='row', metric='cosine', n_jobs=-1):
    """
    Optimized distance matrix calculation with adaptive thresholds and better memory management.
    Specifically optimized for large matrices (20K+ rows).
    """
    try:
        # Ensure data is proper numpy array
        if not isinstance(data, np.ndarray):
            data = np.array(data, dtype=np.float64)
        
        # Get effective number of jobs
        n_jobs = effective_n_jobs(n_jobs)
        
        # Work with the appropriate axis
        work_data = data.T if axis == 'col' else data
        n_samples = work_data.shape[0]
        
        print(f"Computing {axis} distances for {n_samples} samples with {n_jobs} jobs")
        
        # Handle zero-norm vectors for cosine distance
        if metric == 'cosine':
            work_data = _handle_zero_norm_vectors(work_data)
        
        # More aggressive parallel threshold for large matrices
        if n_samples >= 10000:
            # For very large matrices, always use parallel
            distances = _parallel_pdist_large_matrix(work_data, metric, n_jobs)
        elif n_samples >= 2000:
            # For medium-large matrices, use parallel with smaller blocks
            distances = _parallel_pdist(work_data, metric, n_jobs)
        else:
            # Use scipy's optimized pdist for smaller matrices
            distances = pdist(work_data, metric=metric)
        
        # Final validation and cleanup
        distances = _clean_distances(distances)
        
        return distances
        
    except Exception as e:
        print(f"❌ Error in optimized_distance_matrix: {e}")
        traceback.print_exc()
        # Return fallback
        n = data.shape[0] if axis == 'row' else data.shape[1]
        return np.zeros(n * (n - 1) // 2) if n > 1 else np.array([])


def _handle_zero_norm_vectors(data):
    """Handle zero-norm vectors that cause issues with cosine distance."""
    norms = np.linalg.norm(data, axis=1)
    zero_indices = np.where(norms < 1e-12)[0]
    
    if len(zero_indices) > 0:
        print(f"⚠️ Found {len(zero_indices)} zero-norm vectors. Adding small noise...")
        # Create a copy to avoid modifying original data
        data = data.copy()
        noise_scale = 1e-8
        for idx in zero_indices:
            data[idx] = np.random.normal(0, noise_scale, data.shape[1])
    
    return data


def _parallel_pdist_large_matrix(data, metric, n_jobs):
    """
    Specialized parallel pdist for very large matrices (10K+ samples).
    Uses memory-efficient chunking to handle matrices like 20K x 1K.
    """
    n = len(data)
    print(f"Using large matrix optimization for {n} samples")
    
    # For very large matrices, use smaller blocks to manage memory
    # Target ~2MB per block for the data itself
    bytes_per_sample = data.shape[1] * 8  # 8 bytes per float64
    target_block_mb = 2
    target_block_samples = max(100, min(1000, int(target_block_mb * 1024 * 1024 / bytes_per_sample)))
    
    # Create blocks
    n_blocks = max(1, (n + target_block_samples - 1) // target_block_samples)
    block_size = (n + n_blocks - 1) // n_blocks
    
    print(f"Using {n_blocks} blocks of ~{block_size} samples each")
    
    blocks = []
    for i in range(0, n, block_size):
        end_idx = min(i + block_size, n)
        blocks.append(data[i:end_idx])
    
    # Process blocks with memory-conscious batching
    # Limit concurrent jobs to prevent memory explosion
    effective_jobs = min(n_jobs, 4)  # Cap at 4 for memory reasons
    
    def compute_single_block_pair(args):
        i, j = args
        return _compute_block_distances(blocks[i], blocks[j], metric, i == j)
    
    # Generate all block pairs
    block_pairs = [(i, j) for i in range(len(blocks)) for j in range(i, len(blocks))]
    
    # Process in smaller batches to control memory usage
    batch_size = max(1, effective_jobs * 2)
    all_results = []
    
    for batch_start in range(0, len(block_pairs), batch_size):
        batch_end = min(batch_start + batch_size, len(block_pairs))
        batch_pairs = block_pairs[batch_start:batch_end]
        
        batch_results = Parallel(n_jobs=effective_jobs, backend='threading')(
            delayed(compute_single_block_pair)(pair) for pair in batch_pairs
        )
        all_results.extend(batch_results)
        
        # Print progress for very large matrices
        if len(block_pairs) > 20:
            progress = (batch_end / len(block_pairs)) * 100
            print(f"Progress: {progress:.1f}% ({batch_end}/{len(block_pairs)} block pairs)")
    
    # Combine results efficiently
    return _combine_block_results_large(all_results, blocks, n)


def _parallel_pdist(data, metric, n_jobs):
    """Efficient parallel pairwise distance computation for medium-sized matrices."""
    n = len(data)
    
    # Optimize block size for memory and cache efficiency
    # Aim for blocks that fit in L3 cache (~8MB per core)
    target_block_elements = min(1000, max(50, int(np.sqrt(8000000 / (data.shape[1] * 8)))))
    n_blocks = max(1, min(n_jobs * 2, n // target_block_elements))
    block_size = max(1, n // n_blocks)
    
    # Create blocks with proper handling of remainder
    blocks = []
    for i in range(0, n, block_size):
        end_idx = min(i + block_size, n)
        blocks.append(data[i:end_idx])
    
    # Compute block pairs
    block_pairs = [(i, j) for i in range(len(blocks)) for j in range(i, len(blocks))]
    
    # Parallel computation with optimized batch size
    batch_size = max(1, len(block_pairs) // (n_jobs * 4))
    
    def compute_block_batch(batch):
        return [_compute_block_distances(blocks[i], blocks[j], metric, i == j) 
                for i, j in batch]
    
    # Process in batches to reduce memory overhead
    all_results = []
    for batch_start in range(0, len(block_pairs), batch_size):
        batch_end = min(batch_start + batch_size, len(block_pairs))
        batch = block_pairs[batch_start:batch_end]
        
        batch_results = Parallel(n_jobs=n_jobs, backend='threading')(
            delayed(compute_block_batch)([pair]) for pair in batch
        )
        all_results.extend([result[0] for result in batch_results])
    
    # Efficiently combine results into condensed format
    return _combine_block_results(all_results, blocks, n)


def _compute_block_distances(block_i, block_j, metric, is_diagonal):
    """Compute distances between two blocks."""
    try:
        if is_diagonal:
            # For diagonal blocks, use pdist (more efficient)
            return pdist(block_i, metric=metric)
        else:
            # For off-diagonal blocks, compute rectangular distance matrix
            from scipy.spatial.distance import cdist
            return cdist(block_i, block_j, metric=metric)
    except Exception as e:
        print(f"⚠️ Block computation failed: {e}")
        # Return zeros as fallback
        if is_diagonal:
            n = len(block_i)
            return np.zeros(n * (n - 1) // 2)
        else:
            return np.zeros((len(block_i), len(block_j)))


def _combine_block_results(results, blocks, total_n):
    """Efficiently combine block results into condensed distance matrix."""
    condensed_size = total_n * (total_n - 1) // 2
    condensed = np.zeros(condensed_size)
    
    result_idx = 0
    condensed_idx = 0
    
    for i in range(len(blocks)):
        block_i_size = len(blocks[i])
        block_i_start = sum(len(blocks[k]) for k in range(i))
        
        for j in range(i, len(blocks)):
            if i == j:
                # Diagonal block - copy directly to condensed form
                block_result = results[result_idx]
                end_idx = condensed_idx + len(block_result)
                condensed[condensed_idx:end_idx] = block_result
                condensed_idx = end_idx
            else:
                # Off-diagonal block - need to extract relevant parts
                block_result = results[result_idx]
                block_j_start = sum(len(blocks[k]) for k in range(j))
                
                # Convert rectangular block to condensed form contributions
                for ii in range(block_i_size):
                    global_i = block_i_start + ii
                    for jj in range(len(blocks[j])):
                        global_j = block_j_start + jj
                        if global_i < global_j:  # Only upper triangle
                            # Calculate condensed index
                            idx = total_n * global_i - global_i * (global_i + 1) // 2 + global_j - global_i - 1
                            condensed[idx] = block_result[ii, jj]
            
            result_idx += 1
    
    return condensed


# def _combine_block_results_large(results, blocks, total_n):
#     """
#     Memory-efficient combination of block results for very large matrices.
#     Builds condensed matrix incrementally to minimize peak memory usage.
#     """
#     condensed_size = total_n * (total_n - 1) // 2
#     condensed = np.zeros(condensed_size, dtype=np.float32)  # Use float32 to save memory
    
#     print(f"Combining results into condensed matrix of size {condensed_size}")
    
#     result_idx = 0
    
#     # Pre-calculate block start positions
#     block_starts = [0]
#     for i in range(len(blocks)):
#         block_starts.append(block_starts[-1] + len(blocks[i]))
    
#     for i in range(len(blocks)):
#         block_i_start = block_starts[i]
#         block_i_size = len(blocks[i])
        
#         for j in range(i, len(blocks)):
#             block_j_start = block_starts[j]
#             block_j_size = len(blocks[j])
            
#             block_result = results[result_idx]
            
#             if i == j:
#                 # Diagonal block - direct copy to condensed form
#                 # Calculate starting position in condensed array for this block
#                 start_row = block_i_start
#                 start_condensed_idx = start_row * total_n - start_row * (start_row + 1) // 2
                
#                 # Copy the condensed block result
#                 end_idx = start_condensed_idx + len(block_result)
#                 condensed[start_condensed_idx:end_idx] = block_result.astype(np.float32)
                
#             else:
#                 # Off-diagonal block - map rectangular to condensed
#                 for ii in range(block_i_size):
#                     global_i = block_i_start + ii
#                     start_j = max(block_j_start, global_i + 1)  # Only upper triangle
                    
#                     if start_j < block_j_start + block_j_size:
#                         # Calculate condensed indices for this row
#                         local_j_start = start_j - block_j_start
#                         local_j_end = block_j_size
                        
#                         for local_j in range(local_j_start, local_j_end):
#                             global_j = block_j_start + local_j
#                             condensed_idx = total_n * global_i - global_i * (global_i + 1) // 2 + global_j - global_i - 1
#                             condensed[condensed_idx] = block_result[ii, local_j]
            
#             result_idx += 1
    
#     return condensed.astype(np.float64)  # Convert back to float64 for compatibility

def _combine_block_results_large(results, blocks, total_n):
    """
    Memory-efficient combination of block results for very large matrices.
    Optimized for 20K+ samples with vectorized operations.
    """
    condensed_size = total_n * (total_n - 1) // 2
    condensed = np.zeros(condensed_size, dtype=np.float32)
    
    print(f"Combining results into condensed matrix of size {condensed_size}")
    print("Using optimized vectorized combination...")
    
    result_idx = 0
    completed_pairs = 0
    total_pairs = len(blocks) * (len(blocks) + 1) // 2
    
    # Pre-calculate block start positions
    block_starts = np.cumsum([0] + [len(block) for block in blocks])
    
    for i in range(len(blocks)):
        block_i_start = block_starts[i]
        block_i_size = len(blocks[i])
        
        for j in range(i, len(blocks)):
            block_j_start = block_starts[j]
            block_j_size = len(blocks[j])
            
            block_result = results[result_idx]
            
            if i == j:
                # Diagonal block - use vectorized assignment
                if len(block_result) > 0:
                    # Calculate the starting position in condensed array
                    start_idx = _calculate_condensed_start_idx(block_i_start, total_n)
                    end_idx = start_idx + len(block_result)
                    condensed[start_idx:end_idx] = block_result.astype(np.float32)
                
            else:
                # Off-diagonal block - vectorized where possible
                _fill_off_diagonal_block_vectorized(
                    condensed, block_result, 
                    block_i_start, block_i_size,
                    block_j_start, block_j_size, 
                    total_n
                )
            
            result_idx += 1
            completed_pairs += 1
            
            # Progress reporting every 5%
            if completed_pairs % max(1, total_pairs // 20) == 0:
                progress = (completed_pairs / total_pairs) * 100
                print(f"Combination progress: {progress:.1f}% ({completed_pairs}/{total_pairs} block pairs)")
    
    print("✅ Block combination completed")
    return condensed.astype(np.float64)


def _calculate_condensed_start_idx(start_row, total_n):
    """Calculate starting index in condensed array for diagonal block."""
    if start_row == 0:
        return 0
    else:
        return start_row * total_n - start_row * (start_row + 1) // 2


def _fill_off_diagonal_block_vectorized(condensed, block_result, 
                                      block_i_start, block_i_size,
                                      block_j_start, block_j_size, 
                                      total_n):
    """Fill off-diagonal block with vectorized operations where possible."""
    
    # For small blocks, use the old method
    if block_i_size * block_j_size < 10000:
        _fill_off_diagonal_block_simple(condensed, block_result,
                                      block_i_start, block_i_size,
                                      block_j_start, block_j_size, total_n)
        return
    
    # For larger blocks, process row by row with vectorization
    for ii in range(block_i_size):
        global_i = block_i_start + ii
        
        # Only process upper triangle
        start_j = max(block_j_start, global_i + 1)
        
        if start_j < block_j_start + block_j_size:
            local_j_start = start_j - block_j_start
            local_j_end = block_j_size
            
            # Vectorized index calculation for this row
            global_j_array = np.arange(start_j, block_j_start + block_j_size)
            local_j_array = np.arange(local_j_start, local_j_end)
            
            # Calculate condensed indices vectorized
            condensed_indices = (total_n * global_i - 
                               global_i * (global_i + 1) // 2 + 
                               global_j_array - global_i - 1)
            
            # Assign values vectorized
            condensed[condensed_indices] = block_result[ii, local_j_array].astype(np.float32)


def _fill_off_diagonal_block_simple(condensed, block_result,
                                  block_i_start, block_i_size,
                                  block_j_start, block_j_size, total_n):
    """Simple loop-based filling for small blocks."""
    for ii in range(block_i_size):
        global_i = block_i_start + ii
        start_j = max(block_j_start, global_i + 1)
        
        if start_j < block_j_start + block_j_size:
            local_j_start = start_j - block_j_start
            local_j_end = block_j_size
            
            for local_j in range(local_j_start, local_j_end):
                global_j = block_j_start + local_j
                condensed_idx = total_n * global_i - global_i * (global_i + 1) // 2 + global_j - global_i - 1
                condensed[condensed_idx] = block_result[ii, local_j]


def _clean_distances(distances):
    """Clean up distance matrix to ensure validity."""
    if not np.isfinite(distances).all():
        n_bad = np.sum(~np.isfinite(distances))
        print(f"⚠️ Found {n_bad} non-finite distances. Cleaning...")
        distances = np.nan_to_num(distances, nan=1.0, posinf=2.0, neginf=0.0)
    
    # Ensure non-negative distances
    distances = np.maximum(distances, 0.0)
    
    return distances


def optimized_sort_rank_nodes(mat, axis, rank_type):
    """Optimized node ranking using vectorized operations."""
    if axis == 'row':
        data = mat
    else:
        data = mat.T
    
    if rank_type == 'sum':
        ranks = np.sum(data, axis=1)
    elif rank_type == 'var':
        ranks = np.var(data, axis=1)
    else:
        raise ValueError(f"Unknown rank_type: {rank_type}")
    
    # Get sort indices (ascending order)
    sort_indices = np.argsort(ranks)
    
    # Create ranking array where rank[i] is the rank of original element i
    ranking = np.empty_like(sort_indices)
    ranking[sort_indices] = np.arange(len(sort_indices))
    
    return ranking.tolist()


def clust_and_group(net, inst_dm, axis, mat, dist_type='cosine', linkage_type='average',
                    clust_library='scipy', min_samples=1, min_cluster_size=2):
    """Optimized clustering with better error handling."""
    
    try:
        if clust_library == 'scipy':
            Y = hier.linkage(inst_dm, method=linkage_type)
            
        elif clust_library == 'fastcluster':
            import fastcluster
            Y = fastcluster.linkage(inst_dm, method=linkage_type)
            
        elif clust_library == 'hdbscan':
            return _hdbscan_clustering(net, axis, mat, dist_type, 
                                     min_samples, min_cluster_size)
        
        # Generate dendrogram
        Z = hier.dendrogram(Y, no_plot=True)
        inst_clust_order = Z['leaves']
        
        # Generate groups efficiently
        groups = _generate_cluster_groups(Y, inst_dm)
        
        return inst_clust_order, Y, groups
        
    except Exception as e:
        print(f"❌ Clustering failed: {e}")
        traceback.print_exc()
        
        # Return fallback clustering (original order)
        n = len(net.dat['nodes'][axis])
        fallback_order = list(range(n))
        fallback_Y = np.zeros((n-1, 4))  # Minimal linkage matrix
        fallback_groups = {'00': [1] * n}  # All in one group
        
        return fallback_order, fallback_Y, fallback_groups


def _hdbscan_clustering(net, axis, mat, dist_type, min_samples, min_cluster_size):
    """Optimized HDBSCAN clustering."""
    import hdbscan
    from sklearn.decomposition import PCA
    import umap
    
    # Determine data matrix
    if axis == 'row':
        data = mat
    else:
        data = mat.T
    
    # Efficient dimensionality reduction
    n_components = min(50, data.shape[1] - 1, data.shape[0] - 1)
    
    if data.shape[1] > n_components:
        # Use randomized PCA for speed
        pca = PCA(n_components=n_components, random_state=42, svd_solver='randomized')
        low_d_data = pca.fit_transform(data)
    else:
        low_d_data = data
    
    # UMAP with optimized parameters
    umap_reducer = umap.UMAP(
        metric=dist_type,
        n_neighbors=min(15, max(2, len(data) // 10)),
        min_dist=0.0,
        n_components=2,
        random_state=42,
        n_jobs=1  # UMAP threading can conflict with outer parallelization
    )
    
    umap_data = umap_reducer.fit_transform(low_d_data)
    
    # Store UMAP results
    umap_df = pd.DataFrame(
        umap_data.T,
        index=['x', 'y'],
        columns=net.dat['nodes'][axis]
    )
    net.umap[axis] = umap_df
    
    # HDBSCAN clustering
    clusterer = hdbscan.HDBSCAN(
        min_samples=min_samples,
        min_cluster_size=min_cluster_size,
        core_dist_n_jobs=1
    )
    
    clusterer.fit(umap_data)
    
    # Convert to hierarchical format
    Y = clusterer.single_linkage_tree_.to_numpy()
    Z = hier.dendrogram(Y, no_plot=True)
    
    groups = {'hdbscan_clusters': clusterer.labels_.tolist()}
    
    return Z['leaves'], Y, groups


# def _generate_cluster_groups(linkage_matrix, distance_matrix):
#     """Generate cluster groups at different cutoff levels - using original 11 levels."""
#     groups = {}
#     max_dist = distance_matrix.max() if len(distance_matrix) > 0 else 1.0
    
#     # Original cutoff levels (0.0, 0.1, 0.2, ..., 1.0)
#     cutoff_levels = [i / 10.0 for i in range(11)]
    
#     for cutoff in cutoff_levels:
#         try:
#             key = str(int(cutoff * 10)) if cutoff != 1.0 else '10'
#             threshold = cutoff * max_dist
#             groups[key] = hier.fcluster(linkage_matrix, threshold, 'distance').tolist()
#         except Exception as e:
#             print(f"⚠️ Failed to generate groups at cutoff {cutoff}: {e}")
#             # Fallback: single group
#             n_samples = len(linkage_matrix) + 1
#             groups[key] = [1] * n_samples
    
#     return groups

def _generate_cluster_groups(linkage_matrix, distance_matrix):
    """Generate cluster groups at different cutoff levels - using original 11 levels."""
    groups = {}
    max_dist = distance_matrix.max() if len(distance_matrix) > 0 else 1.0
    
    # Original cutoff levels (0.0, 0.1, 0.2, ..., 1.0)
    cutoff_levels = [i / 10.0 for i in range(11)]
    
    for cutoff in cutoff_levels:
        try:
            # Format key to match what viz_json expects: '00', '01', '02', ..., '10'
            # This matches the format created by str(cutoff).replace('.', '')
            key = str(cutoff).replace('.', '')
            if len(key) == 1:  # Handle single digits like '0' -> '00'
                key = '0' + key
            
            threshold = cutoff * max_dist
            groups[key] = hier.fcluster(linkage_matrix, threshold, 'distance').tolist()
        except Exception as e:
            print(f"⚠️ Failed to generate groups at cutoff {cutoff}: {e}")
            # Fallback: single group
            n_samples = len(linkage_matrix) + 1
            groups[key] = [1] * n_samples
    
    return groups




# def group_cutoffs():
#     """Return the original group cutoffs."""
#     return [float(i) / 10 for i in range(11)]

def group_cutoffs():
  all_dist = []
  for i in range(11):
    all_dist.append(float(i) / 10)
  return all_dist


# Legacy functions for backward compatibility
def calc_distance_matrix(tmp_mat, axis, dist_type='cosine'):
    """Legacy function - redirects to optimized version."""
    return optimized_distance_matrix(tmp_mat, axis, dist_type)


def parallel_distance_matrix(data, axis='row', metric='cosine', n_jobs=4):
    """Legacy function - redirects to optimized version."""
    return optimized_distance_matrix(data, axis, metric, n_jobs)


def calculate_blockwise_distances(block_i, block_j, metric='cosine'):
    """Legacy function - redirects to optimized version."""
    return _compute_block_distances(block_i, block_j, metric, block_i is block_j)


def sort_rank_nodes(net, rowcol, rank_type):
    """Legacy function - redirects to optimized version."""
    mat = np.array(net.dat['mat'])
    return optimized_sort_rank_nodes(mat, rowcol, rank_type)