import pandas as pd
import numpy as np
from copy import deepcopy

# def run_norm(net, df=None, norm_type='zscore', axis='row', z_clip=None):
#   '''
#   A dataframe can be passed to run_norm and a normalization will be run (
#   e.g. zscore) on either the rows or columns
#   '''

#   if df is None:
#     df = net.dat_to_df()

#   if norm_type == 'zscore':
#     df, ser_mean, ser_std = zscore_df(df, axis, z_clip=z_clip)

#     net.dat['pre_zscore'] = {}
#     net.dat['pre_zscore']['mean'] = ser_mean.values.tolist()
#     net.dat['pre_zscore']['std'] = ser_std.values.tolist()

#   if norm_type == 'qn':
#     df = qn_df(df, axis)

#   if norm_type == 'umi':
#     df = umi_norm(df)

#   net.df_to_dat(df)

def run_norm(net, df=None, norm_type='zscore', axis='row', z_clip=None):
    '''Debug NaN handling in normalization'''
    
    if df is None:
        df = net.dat_to_df()
        
    
    print(f"🔍 NaN HANDLING DEBUG:")
    print(f"   Input DataFrame shape: {df.shape}")
    
    # Analyze NaN pattern before normalization
    print(f"🔍 NaN analysis BEFORE normalization:")
    total_nans_before = df.isna().sum().sum()
    rows_with_nans = df.isna().any(axis=1).sum()
    cols_with_nans = df.isna().any(axis=0).sum()
    rows_all_nans = df.isna().all(axis=1).sum()
    cols_all_nans = df.isna().all(axis=0).sum()
    
    print(f"   Total NaN values: {total_nans_before}")
    print(f"   Rows with any NaN: {rows_with_nans}")
    print(f"   Cols with any NaN: {cols_with_nans}")
    print(f"   Rows with ALL NaN: {rows_all_nans}")
    print(f"   Cols with ALL NaN: {cols_all_nans}")
    
    if axis == 'row':
        # Check how many valid values each row has for mean/std calculation
        valid_counts_per_row = df.count(axis=1)  # Count non-NaN values per row
        rows_with_insufficient_data = (valid_counts_per_row <= 1).sum()
        print(f"   Rows with ≤1 valid values: {rows_with_insufficient_data}")
        
        if rows_with_insufficient_data > 0:
            insufficient_rows = valid_counts_per_row[valid_counts_per_row <= 1].index.tolist()
            print(f"   Insufficient data row indices: {insufficient_rows[:10]}")
    
    # Perform normalization
    if norm_type == 'zscore':
        print(f"🔍 Performing Z-score normalization...")
        df_normalized, ser_mean, ser_std = zscore_df_with_nan_debug(df, axis, z_clip=z_clip)
        
        net.dat['pre_zscore'] = {}
        net.dat['pre_zscore']['mean'] = ser_mean.values.tolist()
        net.dat['pre_zscore']['std'] = ser_std.values.tolist()
        
    elif norm_type == 'qn':
        df_normalized = qn_df(df, axis)
    elif norm_type == 'umi':
        df_normalized = umi_norm(df)
    else:
        df_normalized = df
    
    # Analyze NaN pattern after normalization
    print(f"🔍 NaN analysis AFTER normalization:")
    total_nans_after = df_normalized.isna().sum().sum()
    rows_with_nans_after = df_normalized.isna().any(axis=1).sum()
    rows_all_nans_after = df_normalized.isna().all(axis=1).sum()
    
    print(f"   Total NaN values: {total_nans_after}")
    print(f"   Rows with any NaN: {rows_with_nans_after}")
    print(f"   Rows with ALL NaN: {rows_all_nans_after}")
    
    # Check if any rows became all NaN after normalization
    newly_all_nan_rows = rows_all_nans_after - rows_all_nans
    if newly_all_nan_rows > 0:
        print(f"   🚨 NEW rows that became ALL NaN: {newly_all_nan_rows}")
    
    # Store pre-conversion state
    shape_before_conversion = df_normalized.shape
    
    # Debug the df_to_dat conversion
    print(f"🔍 Converting to network format...")
    net.df_to_dat(df_normalized)
    
    # Check what happened
    final_shape = np.array(net.dat['mat']).shape
    rows_lost = shape_before_conversion[0] - final_shape[0]
    
    print(f"🔍 CONVERSION RESULTS:")
    print(f"   Before conversion: {shape_before_conversion}")
    print(f"   After conversion: {final_shape}")
    print(f"   Rows lost: {rows_lost}")
    
    if rows_lost > 0:
        print(f"   🚨 {rows_lost} rows were removed during df_to_dat conversion!")
        print(f"   This suggests df_to_dat removes rows with all NaN or invalid values")


def zscore_df_with_nan_debug(df, axis='row', z_clip=None):
    '''Fixed version that handles zero-variance genes properly'''
    
    print(f"🔍 Z-score with NaN handling:")
    
    if axis == 'row':
        df = df.transpose()

    # Check what pandas.mean() and pandas.std() will do with NaN
    print(f"   Calculating statistics (skipna=True by default)...")
    ser_mean = df.mean()  # skipna=True by default
    ser_std = df.std()    # skipna=True by default
    
    # Check for problematic statistics
    nan_means = ser_mean.isna().sum()
    nan_stds = ser_std.isna().sum()
    zero_stds = (ser_std == 0).sum()
    
    print(f"   Columns with NaN mean: {nan_means}")
    print(f"   Columns with NaN std: {nan_stds}")
    print(f"   Columns with zero std: {zero_stds}")
    
    if nan_means > 0:
        nan_mean_cols = ser_mean[ser_mean.isna()].index.tolist()
        print(f"   NaN mean columns: {nan_mean_cols[:5]}")
        
    if nan_stds > 0:
        nan_std_cols = ser_std[ser_std.isna()].index.tolist()
        print(f"   NaN std columns: {nan_std_cols[:5]}")
    
    # 🔧 FIX: Handle zero standard deviation to prevent all-NaN rows
    if zero_stds > 0:
        print(f"   🔧 FIXING {zero_stds} zero-variance columns to prevent NaN creation")
        zero_std_cols = ser_std[ser_std == 0].index.tolist()
        print(f"   Zero std column names: {zero_std_cols[:5]}")
        
        # Create a copy of std series and replace zeros with 1.0
        # This will make the Z-score calculation result in 0.0 for zero-variance genes
        ser_std_fixed = ser_std.copy()
        ser_std_fixed[ser_std == 0] = 1.0
        
        # Perform Z-score calculation with fixed std
        df_z = (df - ser_mean) / ser_std_fixed
        
        # Set zero-variance columns to exactly 0.0 (proper Z-score for no variance)
        zero_var_mask = (ser_std == 0)
        df_z.loc[:, zero_var_mask] = 0.0
        
        print(f"   ✅ Zero-variance columns set to Z-score = 0.0")
        
    else:
        # No zero-variance columns, proceed with normal calculation
        df_z = (df - ser_mean) / ser_std
    
    # Check what happened to NaN values
    nans_before = df.isna().sum().sum()
    nans_after = df_z.isna().sum().sum()
    
    print(f"   NaNs before Z-score: {nans_before}")
    print(f"   NaNs after Z-score: {nans_after}")
    print(f"   Additional NaNs created: {nans_after - nans_before}")
    
    # Check for infinite values created by division by zero
    infs_created = np.isinf(df_z.values).sum()
    print(f"   Infinite values created: {infs_created}")
    
    # Verify no all-NaN rows were created
    all_nan_rows = df_z.isna().all(axis=0).sum()  # axis=0 because we're transposed
    print(f"   All-NaN columns after Z-score: {all_nan_rows}")

    if axis == 'row':
        df_z = df_z.transpose()

    if z_clip is not None:
        df_z = z_clip_fun(df_z, lower=-z_clip, upper=z_clip)

    return df_z, ser_mean, ser_std


# Also update the original zscore_df function for consistency
def zscore_df(df, axis='row', z_clip=None):
    '''Fixed version of the original zscore_df function'''
    
    if axis == 'row':
        df = df.transpose()

    ser_mean = df.mean()
    ser_std = df.std()
    
    # 🔧 FIX: Handle zero standard deviation
    zero_var_mask = (ser_std == 0)
    if zero_var_mask.any():
        print(f"🔧 Handling {zero_var_mask.sum()} zero-variance genes")
        
        # Replace zero std with 1.0 to avoid division by zero
        ser_std_fixed = ser_std.copy()
        ser_std_fixed[zero_var_mask] = 1.0
        
        # Calculate Z-scores
        df_z = (df - ser_mean) / ser_std_fixed
        
        # Set zero-variance genes to Z-score = 0.0
        df_z.loc[:, zero_var_mask] = 0.0
    else:
        # Normal calculation
        df_z = (df - ser_mean) / ser_std

    if axis == 'row':
        df_z = df_z.transpose()

    if z_clip is not None:
        df_z = z_clip_fun(df_z, lower=-z_clip, upper=z_clip)

    return df_z, ser_mean, ser_std

def qn_df(df, axis='row'):
  '''
  do quantile normalization of a dataframe dictionary, does not write to net
  '''
  # using transpose to do row qn
  if axis == 'row':
    df = df.transpose()

  missing_values = df.isnull().values.any()

  # make mask of missing values
  if missing_values:

    # get nan mask
    missing_mask = pd.isnull(df)

    # tmp fill in na with zero, will not affect qn
    df = df.fillna(value=0)

  # calc common distribution
  common_dist = calc_common_dist(df)

  # swap in common distribution
  df = swap_in_common_dist(df, common_dist)

  # swap back in missing values
  if missing_values:
    df = df.mask(missing_mask, other=np.nan)

  # using transpose to do row qn
  if axis == 'row':
    df = df.transpose()

  df_qn = df

  return df_qn

def swap_in_common_dist(df, common_dist):

  col_names = df.columns.tolist()

  qn_arr = np.array([])
  orig_rows = df.index.tolist()

  # loop through each column
  for inst_col in col_names:

    # get the sorted list of row names for the given column
    tmp_series = deepcopy(df[inst_col])
    tmp_series = tmp_series.sort_values(ascending=False)
    sorted_names = tmp_series.index.tolist()

    qn_vect = np.array([])
    for inst_row in orig_rows:
      inst_index = sorted_names.index(inst_row)
      inst_val = common_dist[inst_index]
      qn_vect = np.hstack((qn_vect, inst_val))

    if qn_arr.shape[0] == 0:
      qn_arr = qn_vect
    else:
      qn_arr = np.vstack((qn_arr, qn_vect))

  # transpose (because of vstacking)
  qn_arr = qn_arr.transpose()

  qn_df = pd.DataFrame(data=qn_arr, columns=col_names, index=orig_rows)

  return qn_df

def calc_common_dist(df):
  '''
  calculate a common distribution (for col qn only) that will be used to qn
  '''

  # axis is col
  tmp_arr = np.array([])

  col_names = df.columns.tolist()

  for inst_col in col_names:

    # sort column
    tmp_vect = df[inst_col].sort_values(ascending=False).values

    # stacking rows vertically (will transpose)
    if tmp_arr.shape[0] == 0:
      tmp_arr = tmp_vect
    else:
      tmp_arr = np.vstack((tmp_arr, tmp_vect))

  tmp_arr = tmp_arr.transpose()

  common_dist = tmp_arr.mean(axis=1)

  return common_dist

# def zscore_df(df, axis='row', z_clip=None):
#   '''
#   take the zscore of a dataframe dictionary, does not write to net (self)
#   '''

#   if axis == 'row':
#     df = df.transpose()

#   ser_mean = df.mean()
#   ser_std = df.std()

#   df_z = (df - ser_mean)/ser_std

#   if axis == 'row':
#     df_z = df_z.transpose()

#   if z_clip is not None:
#     df_z = z_clip_fun(df_z, lower=-z_clip, upper=z_clip)

#   return df_z, ser_mean, ser_std

def umi_norm(df):
    # umi norm
    barcode_umi_sum = df.sum()
    df_umi = df.div(barcode_umi_sum)
    return df_umi


def z_clip_fun(df_z, lower=None, upper=None):
  '''
  Trim values at input thresholds using pandas function
  '''
  df_z = df_z.clip(lower=lower, upper=upper)

  return df_z