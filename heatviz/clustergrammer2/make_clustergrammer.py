'''
The clustergrammer python module can be installed using pip:
pip install clustergrammer

or by getting the code from the repo:
https://github.com/MaayanLab/clustergrammer-py
'''

# from clustergrammer import Network
from clustergrammer_fun import Network
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform


def make_cluster(data,assay_type='None',val_scale='None',results_file=False,sampleMetadata=None):

    # if not results_file and assay_type != 'rnaseq':
    #     row_categories = []
    #     if assay_type == 'cytof':
    #         row_categories.append('parent')

    #     if sampleMetadata:
    #         firstKey = list(sampleMetadata.keys())[0]
    #         print(firstKey)
    #         col_categories = list(sampleMetadata[firstKey].keys())
        
    #     column_categories,col_cat_data,row_categories,row_cat_data = csv_to_categories(data,col_categories,row_categories,sampleMetadata)
    

    net = Network()
    net.load_df(data)

    # if val_scale == 'Zscore':
    #     net.normalize(norm_type='zscore', axis='row', keep_orig=False)

    # if not results_file and assay_type != 'rnaseq':
    #     net.add_cats(axis='row', cat_data=row_cat_data)
    #     net.add_cats(axis='col', cat_data=col_cat_data)

    net.cluster(dist_type='cos', dendro=True,
                sim_mat=True, filter_sim=0.1, calc_cat_pval=False, enrichrgram=True)
    
    return net.export_net_json(net_type='viz', indent='no-indent')




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
