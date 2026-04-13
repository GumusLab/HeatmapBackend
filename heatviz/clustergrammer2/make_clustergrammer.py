'''
The clustergrammer python module can be installed using pip:
pip install clustergrammer

or by getting the code from the repo:
https://github.com/MaayanLab/clustergrammer-py
'''

# from clustergrammer import Network
from clustergrammer_fun import Network
import pandas as pd
import cProfile
import pstats
import numpy as np
from scipy.spatial.distance import pdist, squareform


# def make_cluster(data=None,assay_type='None',val_scale='Zscore',results_file=False,sampleMetadata=None):


#     net = Network()
#     # file_content = data.read().decode("utf-8")  # Assuming the file is UTF-8 encoded

#     try:
#         print("Attempting to load file...")
#         net.load_file(data)
#         print("File loaded successfully")
#     except Exception as e:
#         print(f"An error occurred: {e}")
#     if val_scale == 'Zscore':
#         # net.normalize(norm_type='zscore', axis='row', keep_orig=False)
#         net.normalize(norm_type='zscore', axis='row')


#     net.cluster(dist_type='cos', dendro=True,
#                 sim_mat=False, filter_sim=0.1, calc_cat_pval=False, enrichrgram=False)
    
    
#     return net.export_net_json(net_type='viz', indent='no-indent')



def make_cluster(data=None, assay_type='None', val_scale='Zscore', zscore_axis='row', results_file=False,sampleMetadata=None, imputation_method='auto', dist_type='euclidean', linkage_type='average', **imputation_kwargs):

    net = Network()

    try:
        if isinstance(data, pd.DataFrame):
            # NOTE: .copy() is needed - removing it caused issues
            df = data.copy()

            # Create a column list with empty strings for the metadata columns
            # This handles row category columns which appear as empty/NaN in the header
            columns = []
            for i, col in enumerate(df.columns):
                col_str = str(col).lower().strip()
                if 'unnamed' in col_str:
                    # Replace 'Unnamed: X' with empty string (pandas default for unnamed cols)
                    columns.append('')
                elif pd.isna(col) or col == '' or col_str == '' or col_str == 'nan':
                    # Handle NaN, empty strings - these are row category column placeholders
                    columns.append('')
                else:
                    columns.append(col)

            # Convert to a TSV string with blank headers for unnamed columns
            from io import StringIO
            buffer = StringIO()

            # First save as TSV with custom column names
            # Convert to strings for join (handles int column names)
            header_row = '\t'.join(str(c) for c in columns)
            buffer.write(header_row + '\n')

            # Then write the data without header
            data_buffer = StringIO()
            df.to_csv(data_buffer, sep='\t', index=False, header=False)
            buffer.write(data_buffer.getvalue())

            # Load the TSV with proper column handling
            buffer.seek(0)
            net.load_tsv_to_net(buffer)

        elif isinstance(data, str):
            # Assuming 'data' is a string path or raw TSV content
            if "\t" in data or "\n" in data:  # Looks like raw content
                from io import StringIO
                net.load_tsv_to_net(StringIO(data))
            else:
                net.load_file(data)  # Assume it's a file path
        else:
            raise ValueError("Unsupported data type for make_cluster.")

    except Exception as e:
        return {"error": str(e)}

    # Calculate rankings on RAW data (before z-scoring)
    try:
        # Import the ranking function
        from clustergrammer_fun.calc_clust import optimized_sort_rank_nodes
        import numpy as np

        # Get raw matrix as numpy array
        df_raw = net.export_df()
        mat_raw = df_raw.values

        # Calculate rankings for rows (genes/proteins)
        row_rank_sum = optimized_sort_rank_nodes(mat_raw, 'row', 'sum')
        row_rank_var = optimized_sort_rank_nodes(mat_raw, 'row', 'var')

        # Calculate rankings for columns (samples)
        col_rank_sum = optimized_sort_rank_nodes(mat_raw, 'col', 'sum')
        col_rank_var = optimized_sort_rank_nodes(mat_raw, 'col', 'var')

        # Store temporarily on net object for use in clustering
        net.pre_zscore_ranks = {
            'row': {'rank': row_rank_sum, 'rankvar': row_rank_var},
            'col': {'rank': col_rank_sum, 'rankvar': col_rank_var}
        }
    except Exception:
        net.pre_zscore_ranks = None

    if val_scale == 'Zscore' and zscore_axis != 'None':
        net.normalize(norm_type='zscore', axis=zscore_axis)
    
    try:
        # net.cluster(dist_type='cos', dendro=True, sim_mat=False, filter_sim=0.1, clust_library = 'fastcluster', calc_cat_pval=False, enrichrgram=False,imputation_method='knn')
        net.cluster(
            dist_type=dist_type,  # Use the dist_type parameter from function
            linkage_type=linkage_type,  # Use the linkage_type parameter from function
            dendro=True,
            sim_mat=False,
            filter_sim=0.1,
            # clust_library='fastcluster',
            calc_cat_pval=False,
            enrichrgram=False,
            imputation_method=imputation_method,
            **imputation_kwargs
        )
    except Exception as e:
        return {"error": str(e)}

    return net.export_net_json(net_type='viz', indent='no-indent')

# make_cluster()

# net = Network()

# df = pd.read_csv('txt/rnaseq_output.csv', index_col=0)



# # # Check if there are any NaN values in the DataFrame
# # print("Any NaN values in the DataFrame:", df.isnull().values.any())

# # # Option 1: Drop rows with NaN values
# # df_cleaned = df.dropna()

# # # Option 2: Fill NaN values with column mean
# # df_filled = df.fillna(df.mean())

# # # Option 3: Fill NaN values with zero
# # df_filled = df.fillna(0)


# # # Replace infinite values with NaN
# # df_filled.replace([np.inf, -np.inf], np.nan, inplace=True)

# # # Optionally, fill these NaN values with the mean or another strategy
# # df_filled.fillna(df_filled.mean(), inplace=True)



# # load matrix tsv file
# net.load_df(df)

# # # Assuming df is your DataFrame
# # try:
# #     # Compute the pairwise distances
# #     distance_matrix = pdist(df.values, metric='cosine')
# #     print("Max distance:", np.max(distance_matrix))
# #     print("Min distance:", np.min(distance_matrix))
# #     print("Any non-finite values:", not np.isfinite(distance_matrix).all())

# #     # Check for non-finite values
# #     if not np.isfinite(distance_matrix).all():
# #         print("There are non-finite values in the distance matrix.")
# #         indices = np.where(~np.isfinite(distance_matrix))[0]
# #         print("Indices of non-finite values:", indices)

# # except Exception as e:
# #     print("Error while computing distance matrix:", str(e))

# # net.load_file('txt/rc_two_cats.txt')
# # net.load_file('txt/rc_val_cats.txt')

# # optional filtering and normalization
# ##########################################
# # net.filter_sum('row', threshold=20)
# # net.normalize(axis='col', norm_type='zscore', keep_orig=True)
# # net.filter_N_top('row', 250, rank_type='sum')
# # net.filter_threshold('row', threshold=3.0, num_occur=4)
# # net.swap_nan_for_zero()
# # net.downsample(ds_type='kmeans', axis='col', num_samples=10)
# # net.random_sample(random_state=100, num_samples=10, axis='col')
# # net.clip(-6,6)
# # net.filter_cat('row', 1, 'Gene Type: Interesting')
# # net.set_cat_color('col', 1, 'Category: one', 'blue')

# # net.cluster(dist_type='cos',views=['N_row_sum', 'N_row_var'] , dendro=True,
# #              sim_mat=True, filter_sim=0.1, calc_cat_pval=False, enrichrgram=True)

# net.cluster(dist_type='cos', dendro=True,
#                 sim_mat=True, filter_sim=0.1, calc_cat_pval=False, enrichrgram=True)
    

# # write jsons for front-end visualizations
# net.write_json_to_file('viz', 'json/mult_view.json', 'no-indent')
# net.write_json_to_file('sim_row', 'json/mult_view_sim_row.json', 'no-indent')
# net.write_json_to_file('sim_col', 'json/mult_view_sim_col.json', 'no-indent')
