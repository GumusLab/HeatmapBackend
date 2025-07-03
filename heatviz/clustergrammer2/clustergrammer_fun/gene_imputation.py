import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import KNNImputer, IterativeImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
from scipy.stats import pearsonr
from joblib import Parallel, delayed
import warnings
warnings.filterwarnings('ignore')


class GeneExpressionImputer:
    """
    Comprehensive imputation methods for gene expression data.
    
    Methods include:
    - Mean/Median imputation by gene or sample
    - K-Nearest Neighbors (KNN) imputation
    - Iterative imputation (MICE-like)
    - Matrix factorization-based imputation
    - Correlation-based imputation
    - Random Forest imputation
    """
    
    def __init__(self, method='auto', **kwargs):
        """
        Initialize the imputer.
        
        Parameters:
        -----------
        method : str, default='auto'
            Imputation method. Options:
            - 'auto': Automatically select best method based on data characteristics
            - 'mean': Mean imputation by gene
            - 'median': Median imputation by gene
            - 'knn': K-Nearest Neighbors imputation
            - 'iterative': Iterative imputation (MICE-like)
            - 'matrix_factorization': Low-rank matrix completion
            - 'correlation': Correlation-based imputation
            - 'random_forest': Random Forest imputation
            - 'hybrid': Combination of multiple methods
        **kwargs : dict
            Method-specific parameters
        """
        self.method = method
        self.kwargs = kwargs
        self.imputer = None
        self.scaler = None
        self.missing_info = {}
        
    def fit_transform(self, data, gene_names=None, sample_names=None):
        """
        Fit the imputer and transform the data.
        
        Parameters:
        -----------
        data : array-like, shape (n_genes, n_samples)
            Gene expression matrix with potential missing values
        gene_names : list, optional
            Names of genes (rows)
        sample_names : list, optional
            Names of samples (columns)
            
        Returns:
        --------
        imputed_data : np.ndarray
            Imputed gene expression matrix
        """
        data = np.array(data, dtype=np.float64)
        
        # Store original data info
        self.original_shape = data.shape
        self.gene_names = gene_names or [f"Gene_{i}" for i in range(data.shape[0])]
        self.sample_names = sample_names or [f"Sample_{i}" for i in range(data.shape[1])]
        
        # Analyze missing data pattern
        self._analyze_missing_data(data)
        
        # Select imputation method
        if self.method == 'auto':
            self.method = self._select_best_method(data)
            print(f"Auto-selected imputation method: {self.method}")
        
        # Perform imputation
        return self._impute_data(data)
    
    def _analyze_missing_data(self, data):
        """Analyze missing data patterns."""
        missing_mask = np.isnan(data)
        
        self.missing_info = {
            'total_missing': np.sum(missing_mask),
            'missing_percentage': (np.sum(missing_mask) / data.size) * 100,
            'genes_with_missing': np.sum(np.any(missing_mask, axis=1)),
            'samples_with_missing': np.sum(np.any(missing_mask, axis=0)),
            'genes_missing_percentage': (np.sum(missing_mask, axis=1) / data.shape[1]) * 100,
            'samples_missing_percentage': (np.sum(missing_mask, axis=0) / data.shape[0]) * 100
        }
        
        print(f"Missing data analysis:")
        print(f"  Total missing values: {self.missing_info['total_missing']:,}")
        print(f"  Missing percentage: {self.missing_info['missing_percentage']:.2f}%")
        print(f"  Genes with missing values: {self.missing_info['genes_with_missing']}")
        print(f"  Samples with missing values: {self.missing_info['samples_with_missing']}")
    
    def _select_best_method(self, data):
        """Automatically select the best imputation method based on data characteristics."""
        missing_pct = self.missing_info['missing_percentage']
        n_genes, n_samples = data.shape
        
        # Decision tree for method selection
        if missing_pct < 1:
            return 'mean'  # Very few missing values
        elif missing_pct < 5 and n_samples >= 10:
            return 'knn'  # Low missing percentage with enough samples
        elif missing_pct < 15 and n_samples >= 20:
            return 'iterative'  # Moderate missing percentage
        elif missing_pct < 30 and min(n_genes, n_samples) >= 50:
            return 'matrix_factorization'  # High missing percentage with large dataset
        elif n_samples >= 100:
            return 'random_forest'  # Large dataset
        else:
            return 'correlation'  # Fallback for smaller datasets
    
    def _impute_data(self, data):
        """Perform the actual imputation based on selected method."""
        if self.method == 'mean':
            return self._mean_imputation(data)
        elif self.method == 'median':
            return self._median_imputation(data)
        elif self.method == 'knn':
            return self._knn_imputation(data)
        elif self.method == 'iterative':
            return self._iterative_imputation(data)
        elif self.method == 'matrix_factorization':
            return self._matrix_factorization_imputation(data)
        elif self.method == 'correlation':
            return self._correlation_imputation(data)
        elif self.method == 'random_forest':
            return self._random_forest_imputation(data)
        elif self.method == 'hybrid':
            return self._hybrid_imputation(data)
        else:
            raise ValueError(f"Unknown imputation method: {self.method}")
    
    def _mean_imputation(self, data):
        """Mean imputation by gene (row-wise)."""
        print("Performing mean imputation...")
        imputed_data = data.copy()
        
        for i in range(data.shape[0]):
            row = data[i, :]
            if np.any(np.isnan(row)):
                mean_val = np.nanmean(row)
                if np.isnan(mean_val):  # All values are NaN
                    mean_val = 0.0
                imputed_data[i, np.isnan(row)] = mean_val
        
        return imputed_data
    
    def _median_imputation(self, data):
        """Median imputation by gene (row-wise)."""
        print("Performing median imputation...")
        imputed_data = data.copy()
        
        for i in range(data.shape[0]):
            row = data[i, :]
            if np.any(np.isnan(row)):
                median_val = np.nanmedian(row)
                if np.isnan(median_val):  # All values are NaN
                    median_val = 0.0
                imputed_data[i, np.isnan(row)] = median_val
        
        return imputed_data
    
    def _knn_imputation(self, data):
        """K-Nearest Neighbors imputation."""
        print("Performing KNN imputation...")
        
        # KNN works on samples, so we need to transpose for gene expression data
        # where rows are genes and columns are samples
        n_neighbors = self.kwargs.get('n_neighbors', min(5, data.shape[1] - 1))
        n_neighbors = max(1, min(n_neighbors, data.shape[1] - 1))
        
        # Transpose data so samples are rows (required for KNN)
        data_transposed = data.T
        
        imputer = KNNImputer(n_neighbors=n_neighbors, weights='uniform')
        imputed_transposed = imputer.fit_transform(data_transposed)
        
        # Transpose back to original orientation
        return imputed_transposed.T
    
    def _iterative_imputation(self, data):
        """Iterative imputation (MICE-like approach)."""
        print("Performing iterative imputation...")
        
        max_iter = self.kwargs.get('max_iter', 10)
        random_state = self.kwargs.get('random_state', 42)
        
        # Transpose data for iterative imputation
        data_transposed = data.T
        
        imputer = IterativeImputer(
            max_iter=max_iter,
            random_state=random_state,
            initial_strategy='mean',
            n_nearest_features=min(100, data.shape[0])  # Limit features for efficiency
        )
        
        imputed_transposed = imputer.fit_transform(data_transposed)
        return imputed_transposed.T
    
    def _matrix_factorization_imputation(self, data):
        """Low-rank matrix completion using iterative SVD."""
        print("Performing matrix factorization imputation...")
        
        max_iter = self.kwargs.get('max_iter', 50)
        n_components = self.kwargs.get('n_components', min(50, min(data.shape) - 1))
        tolerance = self.kwargs.get('tolerance', 1e-4)
        
        # Initialize with mean imputation
        imputed_data = self._mean_imputation(data)
        missing_mask = np.isnan(data)
        
        for iteration in range(max_iter):
            # Perform SVD
            try:
                U, s, Vt = np.linalg.svd(imputed_data, full_matrices=False)
                
                # Truncate to n_components
                U_trunc = U[:, :n_components]
                s_trunc = s[:n_components]
                Vt_trunc = Vt[:n_components, :]
                
                # Reconstruct matrix
                reconstructed = U_trunc @ np.diag(s_trunc) @ Vt_trunc
                
                # Update only missing values
                old_imputed = imputed_data.copy()
                imputed_data[missing_mask] = reconstructed[missing_mask]
                
                # Check convergence
                if np.allclose(old_imputed[missing_mask], imputed_data[missing_mask], 
                             rtol=tolerance):
                    print(f"Matrix factorization converged after {iteration + 1} iterations")
                    break
                    
            except np.linalg.LinAlgError:
                print("SVD failed, falling back to mean imputation")
                return self._mean_imputation(data)
        
        return imputed_data
    
    def _correlation_imputation(self, data):
        """Correlation-based imputation using gene-gene correlations."""
        print("Performing correlation-based imputation...")
        
        min_correlation = self.kwargs.get('min_correlation', 0.3)
        max_neighbors = self.kwargs.get('max_neighbors', 10)
        
        imputed_data = data.copy()
        missing_mask = np.isnan(data)
        
        # Calculate gene-gene correlation matrix
        corr_matrix = np.corrcoef(data)
        np.fill_diagonal(corr_matrix, 0)  # Remove self-correlation
        
        for gene_idx in range(data.shape[0]):
            if not np.any(missing_mask[gene_idx, :]):
                continue
                
            # Find correlated genes
            correlations = np.abs(corr_matrix[gene_idx, :])
            correlated_genes = np.where(correlations >= min_correlation)[0]
            
            if len(correlated_genes) > max_neighbors:
                # Select top correlated genes
                top_indices = np.argsort(correlations)[-max_neighbors:]
                correlated_genes = top_indices
            
            if len(correlated_genes) == 0:
                # No correlated genes found, use mean
                gene_mean = np.nanmean(data[gene_idx, :])
                if np.isnan(gene_mean):
                    gene_mean = 0.0
                imputed_data[gene_idx, missing_mask[gene_idx, :]] = gene_mean
                continue
            
            # Impute missing values using weighted average of correlated genes
            for sample_idx in np.where(missing_mask[gene_idx, :])[0]:
                weights = []
                values = []
                
                for corr_gene_idx in correlated_genes:
                    if not np.isnan(data[corr_gene_idx, sample_idx]):
                        weights.append(correlations[corr_gene_idx])
                        values.append(data[corr_gene_idx, sample_idx])
                
                if weights:
                    imputed_value = np.average(values, weights=weights)
                    imputed_data[gene_idx, sample_idx] = imputed_value
                else:
                    # No valid correlated values, use gene mean
                    gene_mean = np.nanmean(data[gene_idx, :])
                    if np.isnan(gene_mean):
                        gene_mean = 0.0
                    imputed_data[gene_idx, sample_idx] = gene_mean
        
        return imputed_data
    
    def _random_forest_imputation(self, data):
        """Random Forest-based imputation."""
        print("Performing Random Forest imputation...")
        
        n_estimators = self.kwargs.get('n_estimators', 100)
        max_depth = self.kwargs.get('max_depth', 10)
        random_state = self.kwargs.get('random_state', 42)
        n_jobs = self.kwargs.get('n_jobs', -1)
        
        # Use iterative imputation with Random Forest
        rf_imputer = IterativeImputer(
            estimator=RandomForestRegressor(
                n_estimators=n_estimators,
                max_depth=max_depth,
                random_state=random_state,
                n_jobs=n_jobs
            ),
            max_iter=10,
            random_state=random_state
        )
        
        # Transpose for sample-wise imputation
        data_transposed = data.T
        imputed_transposed = rf_imputer.fit_transform(data_transposed)
        
        return imputed_transposed.T
    
    def _hybrid_imputation(self, data):
        """Hybrid approach combining multiple methods."""
        print("Performing hybrid imputation...")
        
        # Use different methods based on missing data characteristics
        imputed_data = data.copy()
        missing_mask = np.isnan(data)
        
        # For genes with low missing percentage, use correlation-based
        genes_low_missing = np.where(
            (np.sum(missing_mask, axis=1) / data.shape[1]) < 0.1
        )[0]
        
        # For genes with high missing percentage, use matrix factorization
        genes_high_missing = np.where(
            (np.sum(missing_mask, axis=1) / data.shape[1]) >= 0.1
        )[0]
        
        if len(genes_low_missing) > 0:
            # Correlation-based for low missing genes
            corr_imputed = self._correlation_imputation(data[genes_low_missing, :])
            imputed_data[genes_low_missing, :] = corr_imputed
        
        if len(genes_high_missing) > 0:
            # Matrix factorization for high missing genes
            mf_imputed = self._matrix_factorization_imputation(data[genes_high_missing, :])
            imputed_data[genes_high_missing, :] = mf_imputed
        
        return imputed_data
    
    def get_imputation_quality_metrics(self, original_data, imputed_data, 
                                     test_fraction=0.1):
        """
        Evaluate imputation quality by masking known values and comparing predictions.
        
        Parameters:
        -----------
        original_data : np.ndarray
            Original data with missing values
        imputed_data : np.ndarray
            Imputed data
        test_fraction : float
            Fraction of non-missing values to mask for testing
            
        Returns:
        --------
        dict : Quality metrics
        """
        print("Evaluating imputation quality...")
        
        # Find non-missing values in original data
        non_missing_mask = ~np.isnan(original_data)
        non_missing_indices = np.where(non_missing_mask)
        
        if len(non_missing_indices[0]) == 0:
            return {"error": "No non-missing values to evaluate"}
        
        # Randomly select subset for testing
        n_test = int(len(non_missing_indices[0]) * test_fraction)
        test_indices = np.random.choice(len(non_missing_indices[0]), n_test, replace=False)
        
        test_i = non_missing_indices[0][test_indices]
        test_j = non_missing_indices[1][test_indices]
        
        # Mask these values and re-impute
        masked_data = original_data.copy()
        true_values = masked_data[test_i, test_j].copy()
        masked_data[test_i, test_j] = np.nan
        
        # Re-impute with masked data
        temp_imputer = GeneExpressionImputer(method=self.method, **self.kwargs)
        reimputed_data = temp_imputer.fit_transform(masked_data)
        predicted_values = reimputed_data[test_i, test_j]
        
        # Calculate metrics
        mse = np.mean((true_values - predicted_values) ** 2)
        rmse = np.sqrt(mse)
        mae = np.mean(np.abs(true_values - predicted_values))
        
        # Correlation
        correlation, _ = pearsonr(true_values, predicted_values)
        
        # R-squared
        ss_res = np.sum((true_values - predicted_values) ** 2)
        ss_tot = np.sum((true_values - np.mean(true_values)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        return {
            'mse': mse,
            'rmse': rmse,
            'mae': mae,
            'correlation': correlation,
            'r_squared': r_squared,
            'n_test_values': n_test
        }


def update_preprocess_matrix_with_imputation(mat, imputation_method='auto', **kwargs):
    """
    Updated version of _preprocess_matrix function that uses sophisticated imputation.
    
    This function can replace the simple imputation in your existing code.
    """
    # Handle infinite values first
    if not np.isfinite(mat).all():
        print("⚠️ Found infinite values. Clipping to finite range...")
        mat = np.nan_to_num(mat, nan=np.nan, posinf=1e10, neginf=-1e10)
    
    # Check if we have missing values
    if np.isnan(mat).any():
        print(f"⚠️ Found {np.sum(np.isnan(mat))} missing values. Performing imputation...")
        
        # Initialize imputer
        imputer = GeneExpressionImputer(method=imputation_method, **kwargs)
        
        # Perform imputation
        mat = imputer.fit_transform(mat)
        
        # Verify no missing values remain
        if np.isnan(mat).any():
            print("⚠️ Imputation failed to remove all missing values. Using fallback...")
            mat = np.nan_to_num(mat, nan=0.0)
    
    return mat


# Example usage and integration with your existing code
def demo_gene_expression_imputation():
    """Demonstrate various imputation methods on simulated gene expression data."""
    
    # Create synthetic gene expression data
    np.random.seed(42)
    n_genes, n_samples = 1000, 50
    
    # Generate correlated gene expression data
    base_expression = np.random.randn(n_genes, n_samples)
    
    # Add some structure (gene modules)
    for i in range(0, n_genes, 100):
        module_effect = np.random.randn(n_samples) * 2
        base_expression[i:i+100, :] += module_effect
    
    # Introduce missing values (MCAR - Missing Completely At Random)
    missing_rate = 0.1
    missing_mask = np.random.random((n_genes, n_samples)) < missing_rate
    data_with_missing = base_expression.copy()
    data_with_missing[missing_mask] = np.nan
    
    print(f"Created synthetic data: {n_genes} genes, {n_samples} samples")
    print(f"Missing values: {np.sum(missing_mask)} ({missing_rate*100:.1f}%)")
    
    # Test different imputation methods
    methods = ['mean', 'knn', 'iterative', 'matrix_factorization', 'correlation']
    
    results = {}
    for method in methods:
        print(f"\n--- Testing {method} imputation ---")
        
        imputer = GeneExpressionImputer(method=method)
        imputed_data = imputer.fit_transform(data_with_missing)
        
        # Evaluate quality
        quality_metrics = imputer.get_imputation_quality_metrics(
            data_with_missing, imputed_data
        )
        
        results[method] = quality_metrics
        print(f"Quality metrics: {quality_metrics}")
    
    return results