import numpy as np
import time
from scipy.spatial.distance import pdist, cdist, squareform
from joblib import Parallel, delayed
import psutil

def calculate_free_cores():
    # Get the number of logical CPUs
    total_cores = psutil.cpu_count(logical=True)
    
    # Calculate the percentage utilization per CPU core
    cpu_utilization = psutil.cpu_percent(percpu=True)
    
    # Estimate the number of free cores based on low utilization (assuming < 10% as "free")
    free_cores = sum(1 for usage in cpu_utilization if usage < 10)
    
    return free_cores

def calc_distance_matrix(tmp_mat, axis, dist_type='cosine'):
    if axis == 'row':
        inst_dm = pdist(tmp_mat, metric=dist_type)
    elif axis == 'col':
        inst_dm = pdist(tmp_mat.transpose(), metric=dist_type)

    inst_dm[inst_dm < 0] = float(0)

    return inst_dm

def calculate_blockwise_distances(block_i, block_j, metric='cosine'):
    """Calculate distances between rows of block_i and block_j."""
    if block_i is block_j:
        # Diagonal block: Calculate the distance within the block
        return squareform(pdist(block_i, metric=metric))
    else:
        # Off-diagonal block: Calculate the distance between block_i and block_j using cdist
        return cdist(block_i, block_j, metric=metric)

def parallel_distance_matrix(data, metric='cosine', axis='row', n_jobs=-1):
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

# Generate a large test matrix
np.random.seed(42)
large_matrix = np.random.rand(20000, 600)

# Measure the time for the parallelized function
start_time = time.time()
distance_matrix_parallel = parallel_distance_matrix(large_matrix, metric='cosine', axis='row', n_jobs=6)
parallel_time = time.time() - start_time
print(f"Parallel distance matrix calculation took {parallel_time:.4f} seconds.")

# Measure the time for the non-parallelized function
start_time = time.time()
distance_matrix_linear = calc_distance_matrix(large_matrix, axis='row', dist_type='cosine')
linear_time = time.time() - start_time
print(f"Non-parallel distance matrix calculation took {linear_time:.4f} seconds.")

# Check if they match
print("Do the results match?", np.allclose(distance_matrix_parallel, distance_matrix_linear))
