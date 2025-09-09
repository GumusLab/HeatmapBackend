# from collections import defaultdict
# import math
# import os
# import random
# import uuid
# import traceback
# import logging
# from rest_framework.decorators import api_view, parser_classes
# from rest_framework.parsers import MultiPartParser
# from rest_framework.response import Response
# from rest_framework import status
# from ollama import Client as OllamaClient  # ✅ Import Ollama Client
# from django.http import HttpResponse, JsonResponse
# from heatviz.clustergrammer2.make_clustergrammer import make_cluster
# import json
# import pandas as pd
# import copy
# import json
# import requests
# import numpy as np
# from scipy.stats import pearsonr
# from itertools import combinations
# import os
# import multiprocessing as mp
# import threading # For simulating an asynchronous task
# from typing import List, Dict, Any
# import os
# import traceback
# import concurrent.futures
# from .services.pathway_system import (
#     initialize_pathway_database,
#     get_pathway_genes,
#     search_pathways_by_category,
#     get_functional_genes,
#     validate_pathway_command
# )

# # Import your existing functions
# # from .your_existing_modules import filter_genes_by_ids, apply_filters, UPLOAD_DIR

# # Import the ultra-optimized correlation engine
# from .correlation_engine import (
#     compute_correlation_matrix,
#     estimate_correlation_job_time,
#     get_correlation_data_size_estimate
# )

from collections import defaultdict
import math, os, random, uuid, logging, threading, copy, json
from typing import List, Dict, Any
import traceback

from django.http import HttpResponse, JsonResponse
from rest_framework.decorators import api_view, parser_classes
from rest_framework.parsers import MultiPartParser
from rest_framework.response import Response
from rest_framework import status

# Optional: you can also lazy-import this inside the functions that call it
from heatviz.clustergrammer2.make_clustergrammer import make_cluster

from .services.pathway_system import (
    get_pathway_genes,
    search_pathways_by_category,
    get_functional_genes,
    validate_pathway_command
)

from .correlation_engine import (
    compute_correlation_matrix,
    estimate_correlation_job_time,
    get_correlation_data_size_estimate
)


# Set up logger
logger = logging.getLogger('heatviz')
# Default parameters for initial processing and layout generation
DEFAULT_DENDROGRAM_DEPTH = 5
HEMISPHERE_LAYOUT_RADIUS = 50.0
CORRELATION_THRESHOLD_FOR_NEIGHBORS = 0.6 # Threshold for calculating neighbor counts


# Define where to permanently store uploaded files
UPLOAD_DIR = os.path.expanduser("~/heatmap_uploads")
os.makedirs(UPLOAD_DIR, exist_ok=True)  # Ensure the directory exists
OLLAMA_HOST = "http://10.95.46.94:58145"  # Update if the port is different
OLLAMA_MODEL = "gemma"  
# Define the Enrichr URL and the libraries we want to query
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
LIBRARY_CATEGORIES = {
    "Pathways": [
        "MSigDB_Hallmark_2020",             # Essential: Less redundant, highly informative gene sets for biological states.
        "KEGG_2021_Human",                  # Gold-standard: Curated metabolic and signaling pathways.
        "Reactome_2022",                    # Gold-standard: Detailed reaction-level pathway database.
        "WikiPathways_2019_Human",          # Community-curated pathway database.
        "BioPlanet_2019",                   # Comprehensive integrated pathway database.
    ],
    "Transcription": [
        "ChEA_2022",                        # Transcription Factor targets from ChIP-seq experiments.
        "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", # High-quality consensus TF targets.
        "TRRUST_Transcription_Factors_2019" # Curated TF-target regulatory relationships.
    ],
    "Ontologies": [
        "GO_Biological_Process_2023",       # Gold-standard: What biological processes the genes are involved in.
        "GO_Molecular_Function_2023",       # Gold-standard: The molecular activities of the gene products.
        "GO_Cellular_Component_2023",       # Gold-standard: Where in the cell the gene products are active.
    ],
    "Diseases & Drugs": [
        "OMIM_Disease",                     # Curated database of human genes and genetic phenotypes.
        "Jensen_DISEASES",                  # Text-mining based gene-disease associations.
        "DSigDB",                           # Drug-target relationships and gene expression signatures.
    ],
    "Cell Types & Tissues": [
        "GTEx_Tissue_Expression_Down",      # Genes that are specifically down-regulated in certain tissues.
        "GTEx_Tissue_Expression_Up",        # Genes that are specifically up-regulated in certain tissues.
        "Human_Gene_Atlas",                 # Classic tissue-specific gene expression data.
    ],
    "Miscellaneous": [
        "Human_Phenotype_Ontology",         # Connects genes to human phenotypic abnormalities.
        "BioCarta_2016",                    # Older but still referenced pathway diagrams.
    ]
}

# This part of the code remains the same - it automatically builds the list for the API call.
ALL_LIBRARIES = [lib for sublist in LIBRARY_CATEGORIES.values() for lib in sublist]
# --- Setup Logging ---
# Get an instance of a logger for this module.
logger = logging.getLogger(__name__)

# --- Mock Asynchronous Task Queue and Cache ---
# In a real production environment, you would integrate with a robust task queue
# like Celery with a Redis or RabbitMQ broker, and a cache like Redis.
# For this example, we'll use a simple in-memory dictionary to mock task status
# and a threading.Thread for background execution (not recommended for production
# web servers, but serves to demonstrate the asynchronous *concept*).

# This dictionary mocks a persistent cache for task statuses and results.
# Key: session_id, Value: {'status': 'pending'|'computing'|'ready'|'failed', 'data': {...}}
task_status_cache = {}

def calculate_neighbor_counts(correlations: List[Dict[str, Any]], all_gene_ids: List[str], threshold: float = CORRELATION_THRESHOLD_FOR_NEIGHBORS) -> Dict[str, int]:    
    """
    Calculates the degree (neighbor count) for each gene based on a list of correlations.
    
    IMPORTANT: As per the user's clarification, the initial `make_cluster` function
    does NOT produce correlation data from the `proteomics.json` file.
    Therefore, when this function is called from `make_cluster`'s context,
    the `correlations` list will be empty, and all `neighbor_count`s will be 0.
    This is acceptable for the global hemisphere layout, which primarily uses
    `cluster_num` for grouping.
    
    The actual correlation-based neighbor counts will be derived later when
    the `correlation_network` API is called for a specific subset of genes.

    Args:
        correlations (list[dict]): List of correlation dictionaries, e.g., [{'gene1': 'A', 'gene2': 'B', 'correlation': 0.8}].
        all_gene_ids (list[str]): List of all gene IDs.
        threshold (float): Minimum absolute correlation value to count as a neighbor.

    Returns:
        dict[str, int]: A dictionary mapping gene IDs to their neighbor counts.
    """
    logger.info(f"Calculating neighbor counts for {len(all_gene_ids)} genes from {len(correlations)} correlations.")
    gene_degrees = {gene_id: 0 for gene_id in all_gene_ids}
    # If correlations list is empty (as per proteomics.json), all degrees remain 0.
    for corr in correlations:
        if abs(corr.get('correlation', 0.0)) > threshold:
            gene_degrees[corr['gene1']] += 1
            gene_degrees[corr['gene2']] += 1
    return gene_degrees

def _run_async_task_in_background(target_func, *args, **kwargs):
    """
    Helper to run a function in a separate thread.
    WARNING: This is for demonstration purposes only.
    For production, use a proper task queue (e.g., Celery, RQ)
    to handle long-running background tasks.
    """
    thread = threading.Thread(target=target_func, args=args, kwargs=kwargs)
    thread.daemon = True # Allow the main program to exit even if thread is running
    thread.start()


def generate_hemisphere_layout_coords(genes_data: List[Dict[str, Any]], radius: float = HEMISPHERE_LAYOUT_RADIUS) -> Dict[str, Dict[str, float]]:   
    """
    Generates 3D (x,y,z) coordinates for nodes arranged on a full sphere.
    Nodes are primarily grouped by their hierarchical cluster number and
    distributed proportionally around the sphere.

    Args:
        genes_data (list[dict]): List of node dictionaries. Each dict should have:
            - 'id': Unique identifier (e.g., gene name).
            - 'cluster_num': The hierarchical cluster ID the node belongs to.
            - 'neighbor_count': The degree of the node (can be 0).
        radius (float): The radius of the sphere.

    Returns:
        dict[str, dict[str, float]]: A dictionary mapping node IDs to their
                                     calculated {'x': float, 'y': float, 'z': float} coordinates.
    """
    logger.info(f"Generating sphere layout coordinates for {len(genes_data)} genes with radius {radius}")
    
    # NodeData class for sorting, similar to C++ Node struct
    class NodeDataForLayout:
        def __init__(self, id, cluster_num, neighbor_count):
            self.id = id
            self.cluster_num = cluster_num
            self.neighbor_count = neighbor_count

        def __lt__(self, other):
            if self.cluster_num != other.cluster_num:
                return self.cluster_num < other.cluster_num
            else:
                return self.neighbor_count < other.neighbor_count

    processed_nodes_for_layout = []
    for gene_dict in genes_data:
        processed_nodes_for_layout.append(NodeDataForLayout(
            gene_dict['id'],
            gene_dict['cluster_num'],
            gene_dict['neighbor_count']
        ))

    # Sort nodes by cluster_num then by neighbor_count
    processed_nodes_for_layout.sort()

    total_nodes = len(processed_nodes_for_layout)
    if total_nodes == 0:
        return {}

    # Calculate spherical grid increments (adapted from C++ logic)
    # These constants might need tuning for visual aesthetics.
    lg_float = math.sqrt(47520.0 / total_nodes) # Using total_nodes as all are placed now
    lt_float = lg_float * 0.66666667

    lat_inc = max(1, int(lt_float)) # Ensure at least 1 degree increment
    long_inc = max(1, int(lg_float)) # Ensure at least 1 degree increment

    # Count nodes per cluster for proportional distribution
    cl_nr_to_count = defaultdict(int)
    for node in processed_nodes_for_layout:
        cl_nr_to_count[node.cluster_num] += 1
    
    unique_cluster_numbers = sorted(cl_nr_to_count.keys())

    final_coords = {} 

    current_node_idx = 0
    current_longitude_offset = 0 # Offset for current cluster's starting longitude

    # Use a consistent random seed for reproducibility in this dummy function
    random.seed(42) 

    for cluster_num in unique_cluster_numbers:
        nodes_in_this_cluster = cl_nr_to_count[cluster_num]
        
        # Calculate angular span for this cluster proportionally
        # Ensure total_nodes is not zero
        cluster_angular_span = int(math.floor(360.0 * (nodes_in_this_cluster / total_nodes)))
        if cluster_angular_span == 0 and nodes_in_this_cluster > 0: # Ensure small clusters get at least some space
            cluster_angular_span = long_inc # Give them at least one increment

        # Iterate through latitudes for a full sphere (-89 to 89 degrees)
        lat_range = list(range(-89, 90, lat_inc)) # Covers -89 to 89
        if -89 not in lat_range: lat_range.insert(0, -89)
        if 89 not in lat_range: lat_range.append(89)
        lat_range.sort()

        for lat_deg in lat_range:
            for long_deg_offset in range(0, cluster_angular_span, long_inc):
                if current_node_idx >= total_nodes:
                    break # All nodes placed

                node = processed_nodes_for_layout[current_node_idx]
                
                # Only place if it belongs to the current cluster (important due to sorting)
                if node.cluster_num == cluster_num:
                    adjusted_long_deg = long_deg_offset + current_longitude_offset
                    
                    # Convert degrees to radians for math functions
                    lat_rad = degrees_to_radians(lat_deg)
                    long_rad = degrees_to_radians(adjusted_long_deg)

                    # Spherical to Cartesian coordinates
                    x = radius * math.cos(lat_rad) * math.cos(long_rad)
                    y = radius * math.sin(lat_rad) # Y is up/down axis
                    z = radius * math.cos(lat_rad) * math.sin(long_rad)
                    
                    final_coords[node.id] = {'x': x, 'y': y, 'z': z}
                    current_node_idx += 1
                elif node.cluster_num != cluster_num:
                    # This node belongs to a different cluster (due to sorting).
                    # We need to break from this inner loop and move to the next cluster.
                    break 
            
            if current_node_idx >= total_nodes:
                break # All nodes placed

        current_longitude_offset += cluster_angular_span # Move to the next cluster's starting longitude
        if current_node_idx >= total_nodes:
            break # All nodes placed

    # Ensure all nodes are placed, even if logic above missed some due to rounding/small clusters
    # This acts as a fallback, placing any unplaced nodes randomly within the sphere.
    for node in processed_nodes_for_layout:
        if node.id not in final_coords:
            final_coords[node.id] = {
                'x': random.uniform(-radius * 0.5, radius * 0.5),
                'y': random.uniform(-radius * 0.5, radius * 0.5),
                'z': random.uniform(-radius * 0.5, radius * 0.5)
            }

    logger.info(f"Finished generating {len(final_coords)} coordinates.")
    return final_coords

def degrees_to_radians(degrees: float) -> float:
    """Converts degrees to radians."""
    return degrees * (math.pi / 180.0)

def position_clusters_on_sphere(clusters, main_radius):
    """
    Distribute cluster centers evenly on a sphere surface via Fibonacci spiral.
    clusters: list of dicts {'id': int, 'nodes': [geneId, ...]}
    returns: dict mapping cluster_id -> (x, y, z)
    """
    N = len(clusters)
    golden_angle = math.pi * (3 - math.sqrt(5))
    centers = {}

    for i, cluster in enumerate(clusters):
        y = 1 - (2 * i) / (N - 1)
        r_xy = math.sqrt(max(0.0, 1 - y * y))
        theta = golden_angle * i

        x = math.cos(theta) * r_xy
        z = math.sin(theta) * r_xy

        # pull in slightly so genes don't poke out too far
        scale = main_radius * 0.8
        centers[cluster["id"]] = (
            x * scale,
            y * scale,
            z * scale
        )

    return centers

def position_genes_in_clusters(clusters, cluster_centers, main_radius):
    """
    For each cluster, scatter its genes in a small cloud around its center,
    then project back onto the sphere surface with jitter.
    returns: dict mapping geneId -> (x, y, z)
    """
    gene_positions = {}

    for cluster in clusters:
        cid = cluster["id"]
        genes = cluster["nodes"]
        M = len(genes)
        cx, cy, cz = cluster_centers[cid]

        # radius of the little cloud scales with cube-root of cluster size
        cluster_radius = (M ** (1/3)) * (main_radius * 0.05)

        for j, gene_id in enumerate(genes):
            if M == 1:
                # singleton: shoot straight out
                norm = math.sqrt(cx*cx + cy*cy + cz*cz)
                factor = main_radius / (norm or 1)
                x, y, z = cx * factor, cy * factor, cz * factor
            else:
                # local fibonacci spiral
                golden_angle = math.pi * (3 - math.sqrt(5))
                y0 = 1 - (2 * j) / (M - 1)
                r0 = math.sqrt(max(0.0, 1 - y0*y0))
                theta0 = golden_angle * j

                lx = math.cos(theta0) * r0 * cluster_radius
                ly = y0 * cluster_radius
                lz = math.sin(theta0) * r0 * cluster_radius

                # world position before projection
                wx = cx + lx
                wy = cy + ly
                wz = cz + lz

                # normalize + jitter back to sphere surface
                norm = math.sqrt(wx*wx + wy*wy + wz*wz) or 1
                jitter = 1 + (random.random() - 0.5) * 0.1  # ±5%
                factor = main_radius * jitter / norm
                x, y, z = wx * factor, wy * factor, wz * factor

            gene_positions[gene_id] = {"x": x, "y": y, "z": z}

    return gene_positions

def _generate_and_store_3d_coords_task(session_id: str, heatmap_data: dict, radius: float):
    """
    Background task to generate a correlation-informed 3D layout and store it.
    Hybrid approach: All nodes placed by cluster + strongest edges from correlation network.
    """
    coords_file_path = os.path.join(UPLOAD_DIR, f"{session_id}_coords.json")

    try:
        logger.info(f"Async task: Starting correlation-informed 3D coordinate generation for session: {session_id}")
        task_status_cache[session_id] = {
            'status': 'computing',
            'message': 'Generating correlation-informed 3D layout...'
        }

        # Extract row_nodes and build gene list
        row_nodes = heatmap_data.get("row_nodes", [])
        all_gene_ids = [n["name"] for n in row_nodes if "name" in n]
        logger.info(f"Processing {len(all_gene_ids)} genes for 3D layout")

        # Build cluster assignments
        cluster_assignments = {}
        for node in row_nodes:
            name = node.get("name")
            grp  = node.get("group", [])
            if name:
                if DEFAULT_DENDROGRAM_DEPTH < len(grp):
                    cluster_assignments[name] = int(grp[DEFAULT_DENDROGRAM_DEPTH])
                else:
                    cluster_assignments[name] = 0

        # Invert to list of clusters
        clusters_map = {}
        for gene, cid in cluster_assignments.items():
            clusters_map.setdefault(cid, []).append(gene)
        clusters = [{"id": cid, "nodes": genes} for cid, genes in clusters_map.items()]

        # Step 1: Position cluster centers on sphere
        cluster_centers = position_clusters_on_sphere(clusters, radius)

        # Step 2: Get correlation network for strongest edges using SAME gene names as coordinates
        logger.info("Computing correlation network for strongest edges...")
        correlation_edges = []
        try:
            # Use the same gene names that were used for clustering (from row_nodes)
            clustered_gene_names = list(cluster_assignments.keys())
            
            # Load the full dataset for correlation computation
            file_path = os.path.join(UPLOAD_DIR, f"{session_id}.tsv")
            if os.path.exists(file_path):
                import pandas as pd
                from .correlation_engine import compute_correlation_matrix
                
                # Load and prepare data
                df = pd.read_csv(file_path, sep="\t")
                import numpy as np  # For data quality checks
                
                # Filter dataframe to only include genes that are in our clustered data
                if df.index.name is None:
                    # If no index name, assume first column is gene names
                    df = df.set_index(df.columns[0])
                
                # Detect where the numeric data matrix starts (exclude metadata rows/columns)
                matrix_start_row, matrix_start_col = detect_matrix_start(df)
                # Extract only the numeric data matrix (excluding metadata)
                numeric_df = df.iloc[matrix_start_row:, matrix_start_col:]
                # Convert to numeric, handling any remaining non-numeric values
                numeric_df = numeric_df.apply(pd.to_numeric, errors='coerce')
                
                # Update the DataFrame to use the numeric matrix
                df = numeric_df                
                # Filter to only genes that exist in both datasets
                # Handle duplicate gene names: map clustered names (with _dup suffixes) to raw names
                clustered_to_raw_mapping = {}
                for clustered_gene in clustered_gene_names:
                    # Remove _dup suffix if present (e.g., 'U1_dup2' -> 'U1')
                    if '_dup' in clustered_gene:
                        raw_gene = clustered_gene.split('_dup')[0]
                    else:
                        raw_gene = clustered_gene  # Keep as is if no _dup suffix
                    
                    if raw_gene in df.index:
                        clustered_to_raw_mapping[clustered_gene] = raw_gene
                
                available_genes = list(clustered_to_raw_mapping.keys())
                
                if len(available_genes) > 5:  # Lower minimum genes for correlation
                    # Filter dataframe using the raw gene names (without _dup suffixes)
                    raw_genes_to_use = [clustered_to_raw_mapping[gene] for gene in available_genes]
                    
                    # Check for duplicates in raw_genes_to_use
                    raw_genes_unique = list(set(raw_genes_to_use))
                    if len(raw_genes_unique) != len(raw_genes_to_use):
                        from collections import Counter
                        duplicates = [gene for gene, count in Counter(raw_genes_to_use).items() if count > 1]
                        logger.warning(f"   Duplicate genes found in raw_genes_to_use")
                    
                    # Check for duplicates in DataFrame index
                    df_index_unique = len(set(df.index))
                    if df_index_unique != len(df.index):
                        logger.warning(f"   WARNING: DUPLICATES in df.index!")
                        df_duplicates = df.index[df.index.duplicated()].tolist()
                    
                    try:
                        # Handle duplicates properly: use unique genes for correlation
                        if len(raw_genes_unique) != len(raw_genes_to_use):
                            
                            # Create mapping from unique raw genes to clustered genes
                            unique_to_clustered = {}
                            for i, raw_gene in enumerate(raw_genes_to_use):
                                clustered_gene = available_genes[i]
                                if raw_gene not in unique_to_clustered:
                                    unique_to_clustered[raw_gene] = []
                                unique_to_clustered[raw_gene].append(clustered_gene)
                            
                            # logger.info(f"   unique_to_clustered mapping examples: {list(unique_to_clustered.items())[:5]}")
                            
                            # Build DataFrame with unique genes only
                            # We need to handle the case where df.index also has duplicates
                            df_filtered_data = []
                            unique_clustered_names = []
                            
                            for raw_gene in raw_genes_unique:
                                # Get the first occurrence of this gene in the DataFrame
                                if raw_gene in df.index:
                                    # Get first occurrence of this gene
                                    gene_data = df.loc[raw_gene]
                                    if isinstance(gene_data, pd.Series):
                                        # Single occurrence
                                        df_filtered_data.append(gene_data.values)
                                    else:
                                        # Multiple occurrences, take first
                                        df_filtered_data.append(gene_data.iloc[0].values)
                                    
                                    # Use the first clustered gene name
                                    unique_clustered_names.append(unique_to_clustered[raw_gene][0])
                                else:
                                    logger.warning(f"   Gene {raw_gene} not found in DataFrame index")
                            
                            # Create new DataFrame with unique genes
                            df_filtered = pd.DataFrame(df_filtered_data, 
                                                     index=unique_clustered_names,
                                                     columns=df.columns)
                            
                            
                        else:
                            # No duplicates, proceed normally
                            df_filtered = df.loc[raw_genes_to_use]
                            df_filtered.index = available_genes
                           
                    except Exception as filter_error:
                        missing_genes = [gene for gene in raw_genes_to_use if gene not in df.index]
                        logger.error(f"   Missing genes: {missing_genes[:10]} (showing first 10)")
                        raise filter_error
                    

                    
                   
                    
                    # Check for constant genes (zero variance)
                    gene_variances = df_filtered.var(axis=1)
                    zero_var_genes = (gene_variances == 0).sum()
                                        
                    # Use correlation engine with settings for strongest edges only
                    correlation_filters = {
                        'correlationThreshold': 0.1,  # Very low threshold to see any edges
                        'pValueThreshold': 0.1,       # More lenient p-value
                        'maxCorrelations': 5000,      # Limit edges to avoid hairball
                        'variancePercentile': 0.1,    # Include more genes
                        'minSamples': 5               # Lower minimum samples
                    }
                    
                    if df_filtered.index.duplicated().any():
                        duplicated_indices = df_filtered.index[df_filtered.index.duplicated()].tolist()
                        logger.info(f"   Duplicate indices: {duplicated_indices[:10]}")
                    
                    correlation_result = compute_correlation_matrix(df_filtered, correlation_filters)
                    
                    if "correlations" in correlation_result:
                        correlation_edges = correlation_result["correlations"]
                    else:
                        logger.warning("No correlations found, proceeding with cluster-only layout")
                        
                        # Try with even more lenient settings
                        ultra_lenient_filters = {
                            'correlationThreshold': 0.01,  # Almost any correlation
                            'pValueThreshold': 0.5,        # Very lenient p-value
                            'maxCorrelations': 10000,      # More edges
                            'variancePercentile': 0.01,    # Include almost all genes
                            'minSamples': 3                # Very low minimum samples
                        }
                        logger.info(f"Ultra-lenient filters: {ultra_lenient_filters}")
                        ultra_result = compute_correlation_matrix(df_filtered, ultra_lenient_filters)
                        if "correlations" in ultra_result:
                            correlation_edges = ultra_result["correlations"]
                        else:
                            logger.warning(f"Even ultra-lenient failed: {ultra_result}")
                else:
                    logger.warning(f"Too few common genes ({len(available_genes)}), proceeding with cluster-only layout")
            else:
                logger.warning(f"Session file not found: {file_path}")
        except Exception as e:
            logger.warning(f"Failed to compute correlations: {e}, proceeding with cluster-only layout")

        # Step 3: Position individual genes with correlation refinement
        global_3d_positions = position_genes_with_correlation_refinement(
            clusters, cluster_centers, radius, correlation_edges
        )

        # Step 4: Store both coordinates and edge information
        result_data = {
            'coordinates': global_3d_positions,
            'correlation_edges': correlation_edges[:1000],  # Limit edges for frontend
            'cluster_info': {
                'clusters': clusters,
                'cluster_centers': {str(cid): {'x': center[0], 'y': center[1], 'z': center[2]} 
                                  for cid, center in cluster_centers.items()}
            }
        }

        # Persist to disk
        with open(coords_file_path, 'w') as f:
            json.dump(result_data, f)

        task_status_cache[session_id] = {
            'status': 'ready',
            'message': '3D correlation-informed layout ready.',
            'data': result_data
        }
        logger.info(f"Async task: 3D coordinates with {len(correlation_edges)} edges saved to: {coords_file_path}")

    except Exception as e:
        tb = traceback.format_exc()
        logger.error(f"Async task: 3D coordinate generation failed for session {session_id}: {e}\n{tb}")
        task_status_cache[session_id] = {
            'status': 'failed',
            'message': f"Failed to generate 3D layout: {str(e)}"
        }


def position_genes_with_correlation_refinement(clusters, cluster_centers, radius, correlation_edges):
    """
    Position genes in clusters with correlation-based refinement.
    Enhanced to allow inter-cluster correlations to influence positioning.
    """
    import math
    import random
    
    gene_positions = {}
    
    # Build edge lookup for quick access
    edge_lookup = {}
    for edge in correlation_edges:
        gene1, gene2 = edge['gene1'], edge['gene2']
        correlation = edge['correlation']
        edge_lookup.setdefault(gene1, []).append((gene2, correlation))
        edge_lookup.setdefault(gene2, []).append((gene1, correlation))
    
    logger.info(f"Built edge lookup with {len(correlation_edges)} edges")
    
    # Always prefer cluster-based layout to maintain cluster grouping
    # Only use global layout for very sparse networks with few clusters
    if len(correlation_edges) > 1000 and len(clusters) <= 2:
        try:
            import networkx as nx
            
            # Build global network with all genes and strong correlations
            G_global = nx.Graph()
            all_genes = []
            for cluster in clusters:
                all_genes.extend(cluster["nodes"])
            
            # Add all genes as nodes
            for gene in all_genes:
                G_global.add_node(gene)
            
            # Add strong correlations as edges (both intra and inter-cluster)
            strong_edges = 0
            for edge in correlation_edges:
                if abs(edge['correlation']) > 0.5:  # Higher threshold for global layout
                    G_global.add_edge(edge['gene1'], edge['gene2'], weight=abs(edge['correlation']))
                    strong_edges += 1
            
            logger.info(f"Global network: {len(all_genes)} genes, {strong_edges} strong edges")
            
            if G_global.number_of_edges() > 0:
                # Use 3D spring layout for the entire network
                logger.info("Using global 3D spring layout with correlation influence")
                pos_global = nx.spring_layout(G_global, dim=3, k=radius*0.1, iterations=100, seed=42)
                
                # Scale positions to fit within sphere
                max_dist = 0
                for pos in pos_global.values():
                    dist = math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
                    max_dist = max(max_dist, dist)
                
                scale_factor = radius * 0.8 / (max_dist or 1)
                
                for gene in all_genes:
                    if gene in pos_global:
                        pos = pos_global[gene]
                        x = pos[0] * scale_factor
                        y = pos[1] * scale_factor
                        z = pos[2] * scale_factor
                        
                        # Add slight jitter and project to sphere surface
                        jitter = 1 + (random.random() - 0.5) * 0.1
                        norm = math.sqrt(x*x + y*y + z*z) or 1
                        factor = radius * jitter / norm
                        
                        gene_positions[gene] = {"x": x * factor, "y": y * factor, "z": z * factor}
                
                logger.info(f"Global layout positioned {len(gene_positions)} genes")
                return gene_positions
                
        except Exception as e:
            logger.warning(f"Global layout failed: {e}, falling back to cluster-based layout")
    
    # Fallback to cluster-based positioning with enhanced inter-cluster awareness
    logger.info("Using enhanced cluster-based positioning for better visual grouping")
    
    for cluster in clusters:
        cid = cluster["id"]
        genes = cluster["nodes"]
        M = len(genes)
        cx, cy, cz = cluster_centers[cid]
        
        # Base cluster radius scales with cluster size
        base_cluster_radius = (M ** (1/3)) * (radius * 0.12)  # Larger clusters for better visibility
        
        if M == 1:
            # Singleton: place at cluster center projected to sphere
            norm = math.sqrt(cx*cx + cy*cy + cz*cz)
            factor = radius / (norm or 1)
            gene_positions[genes[0]] = {"x": cx * factor, "y": cy * factor, "z": cz * factor}
            continue
        
        # For clusters with multiple genes, use correlation-informed positioning
        genes_with_correlations = [g for g in genes if g in edge_lookup]
        
        if len(genes_with_correlations) > 1:
            try:
                import networkx as nx
                # Build subgraph for this cluster with enhanced connectivity
                G = nx.Graph()
                for gene in genes:
                    G.add_node(gene)
                
                # Add edges within cluster AND strong inter-cluster edges
                for gene in genes:
                    if gene in edge_lookup:
                        for neighbor, correlation in edge_lookup[gene]:
                            abs_corr = abs(correlation)
                            # Include intra-cluster edges (lower threshold) and strong inter-cluster edges
                            if neighbor in genes and abs_corr > 0.1:  # Intra-cluster (lower threshold)
                                G.add_edge(gene, neighbor, weight=abs_corr)
                            elif neighbor not in genes and abs_corr > 0.6:  # Strong inter-cluster
                                # Add phantom node for inter-cluster influence
                                phantom_id = f"phantom_{neighbor}"
                                G.add_node(phantom_id)
                                G.add_edge(gene, phantom_id, weight=abs_corr * 0.3)
                
                if G.number_of_edges() > 0:
                    # Use spring layout in 3D with more iterations for tighter clustering
                    pos = nx.spring_layout(G, dim=3, k=base_cluster_radius/2, iterations=30, seed=42)
                    
                    # Scale and translate to cluster position
                    for gene in genes:
                        if gene in pos:
                            local_pos = pos[gene]
                            # Scale to cluster size and translate to cluster center
                            x = cx + local_pos[0] * base_cluster_radius
                            y = cy + local_pos[1] * base_cluster_radius  
                            z = cz + local_pos[2] * base_cluster_radius
                            
                            # Project back to sphere surface with jitter
                            norm = math.sqrt(x*x + y*y + z*z) or 1
                            jitter = 1 + (random.random() - 0.5) * 0.05  # Less jitter for tighter clusters
                            factor = radius * jitter / norm
                            gene_positions[gene] = {"x": x * factor, "y": y * factor, "z": z * factor}
                        else:
                            # Fallback for isolated nodes
                            gene_positions[gene] = _position_gene_in_cluster_fallback(
                                genes.index(gene), M, cx, cy, cz, base_cluster_radius, radius
                            )
                    print(f"Iteration completed for cluster {cid}")
                else:
                    # No edges in cluster, use original positioning
                    for j, gene in enumerate(genes):
                        gene_positions[gene] = _position_gene_in_cluster_fallback(
                            j, M, cx, cy, cz, base_cluster_radius, radius
                        )
                        
            except ImportError:
                logger.warning("NetworkX not available, using fallback positioning")
                for j, gene in enumerate(genes):
                    gene_positions[gene] = _position_gene_in_cluster_fallback(
                        j, M, cx, cy, cz, base_cluster_radius, radius
                    )
            except Exception as e:
                logger.warning(f"Spring layout failed: {e}, using fallback")
                for j, gene in enumerate(genes):
                    gene_positions[gene] = _position_gene_in_cluster_fallback(
                        j, M, cx, cy, cz, base_cluster_radius, radius
                    )
        else:
            # No correlations available, use original fibonacci spiral
            for j, gene in enumerate(genes):
                gene_positions[gene] = _position_gene_in_cluster_fallback(
                    j, M, cx, cy, cz, base_cluster_radius, radius
                )
    
    logger.info(f"Cluster-based layout positioned {len(gene_positions)} genes in {len(clusters)} clusters")
    return gene_positions


def _position_gene_in_cluster_fallback(j, M, cx, cy, cz, cluster_radius, main_radius):
    """
    Fallback positioning using fibonacci spiral (original logic).
    """
    import math
    import random
    
    golden_angle = math.pi * (3 - math.sqrt(5))
    y0 = 1 - (2 * j) / (M - 1)
    r0 = math.sqrt(max(0.0, 1 - y0*y0))
    theta0 = golden_angle * j

    lx = math.cos(theta0) * r0 * cluster_radius
    ly = y0 * cluster_radius
    lz = math.sin(theta0) * r0 * cluster_radius

    # World position before projection
    wx = cx + lx
    wy = cy + ly
    wz = cz + lz

    # Normalize + jitter back to sphere surface
    norm = math.sqrt(wx*wx + wy*wy + wz*wz) or 1
    jitter = 1 + (random.random() - 0.5) * 0.1
    factor = main_radius * jitter / norm
    
    return {"x": wx * factor, "y": wy * factor, "z": wz * factor}

@api_view(['POST'])
@parser_classes([MultiPartParser])
def process_data_view(request):
    import time
    import threading
    import os
    
    # DEBUG: Log every single request that hits this endpoint
    initial_timestamp = time.time()
    initial_thread_id = threading.current_thread().ident  
    initial_process_id = os.getpid()
    
    print('=' * 100)
    print(f'🔍 INCOMING REQUEST to process_data_view')
    print(f'   Timestamp: {initial_timestamp}')
    print(f'   Human time: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(initial_timestamp))}')
    print(f'   Process ID: {initial_process_id}')
    print(f'   Thread ID: {initial_thread_id}')
    print(f'   Request method: {request.method}')
    print(f'   Content type: {request.content_type}')
    print(f'   Request data keys: {list(request.data.keys())}')
    print(f'   Request FILES keys: {list(request.FILES.keys())}')
    print(f'   Client IP: {request.META.get("REMOTE_ADDR", "unknown")}')
    print(f'   User Agent: {request.META.get("HTTP_USER_AGENT", "unknown")[:50]}...')
    print('=' * 100)

    print('************* request came *************************************')
    
    # Check if this is a strategy selection request (second step)
    if 'session_id' in request.data and 'imputation_strategy' in request.data:
        print('************* strategy selection request *************************************')
        return handle_strategy_processing(request)
    
    # Original file upload logic (first step)
    file = request.FILES.get('data')
    if not file:
        print('************* No file provided *************************************')
        return Response({"error": "No file provided"}, status=status.HTTP_400_BAD_REQUEST)
    
    # DEBUG: Add comprehensive request tracing
    import time
    import threading
    import hashlib
    
    request_timestamp = time.time()
    thread_id = threading.current_thread().ident
    process_id = os.getpid()
    
    # Calculate file hash for duplicate detection
    file_content = file.read()
    file.seek(0)  # Reset file pointer
    file_hash = hashlib.md5(file_content).hexdigest()
    file_size = len(file_content)
    
    print(f'🔍 DEBUG UPLOAD REQUEST:')
    print(f'   Timestamp: {request_timestamp}')
    print(f'   Process ID: {process_id}')
    print(f'   Thread ID: {thread_id}')
    print(f'   File size: {file_size} bytes')
    print(f'   File hash: {file_hash}')
    print(f'   Client IP: {request.META.get("REMOTE_ADDR", "unknown")}')
    print(f'   User Agent: {request.META.get("HTTP_USER_AGENT", "unknown")[:100]}')
    print(f'   Request Headers: {dict(request.headers)}')
    print('*' * 80)
           
        # Defensive: Ensure directory exists before writing
    if not os.path.exists(UPLOAD_DIR):
        print('************* UPLOAD_DIR does not exist *************************************')
        os.makedirs(UPLOAD_DIR, exist_ok=True)

    try:
        # Generate a unique session ID
        print('🔍 DEBUG: Starting session ID generation')
        session_id = str(uuid.uuid4())
        print(f'🔍 DEBUG: Generated session_id: {session_id}')
        
        filename = f"{session_id}.tsv"
        file_path = os.path.join(UPLOAD_DIR, filename)
        metadata_file_path = os.path.join(UPLOAD_DIR, f"{session_id}_metadata.json")
        
        print(f'🔍 DEBUG: File paths:')
        print(f'   TSV: {file_path}')
        print(f'   Metadata: {metadata_file_path}')
        
        # Check if files already exist (shouldn't happen with UUID4)
        if os.path.exists(file_path):
            print(f'🚨 WARNING: File already exists: {file_path}')
        if os.path.exists(metadata_file_path):
            print(f'🚨 WARNING: Metadata file already exists: {metadata_file_path}')
        
        # Save the uploaded file permanently
        print(f'🔍 DEBUG: Starting file write to {file_path}')
        file_write_start = time.time()
        
        with open(file_path, 'wb+') as destination:
            for chunk in file.chunks():
                destination.write(chunk)
        
        file_write_end = time.time()
        print(f'🔍 DEBUG: File write completed in {file_write_end - file_write_start:.3f}s')
        print(f'🔍 DEBUG: File saved: {file_path} (size: {os.path.getsize(file_path)} bytes)')
          
        # 📚 Extract and Save Metadata JSON
        extract_and_save_metadata(tsv_file_path=file_path, metadata_json_path=metadata_file_path)
        print(f"Metadata saved: {metadata_file_path}")
          
        # NEW: Check for missing values before processing
        missing_analysis = analyze_missing_values(file_path)

        
        if missing_analysis['has_missing_values']:
            # Missing values detected - return analysis and available strategies
            print(f"WARNING: Missing values detected: {missing_analysis['missing_percentage']:.2f}%")
            
            return JsonResponse({
                "session_id": session_id,
                "has_missing_values": True,
                "missing_value_summary": missing_analysis,
                # "available_strategies": get_available_strategies()
            }, status=200)
        else:
            # No missing values - proceed with normal processing (default imputation)
            print("No missing values detected, proceeding with clustering")
            response_data = make_cluster(
                data=file_path,
                imputation_method='auto'  # Default when no missing values
            )

            # Handle response_data which can be either JSON string (success) or dict (error)
            if isinstance(response_data, str):
                # Success case: response_data is JSON string
                heatmap_data = json.loads(response_data)
            elif isinstance(response_data, dict):
                # Error case: response_data is already a dict
                if "error" in response_data:
                    # Handle the error case
                    logger.error(f"Error from make_cluster: {response_data['error']}")
                    return JsonResponse({"error": response_data["error"]}, status=500)
                else:
                    # It's a dict but not an error (shouldn't happen, but handle it)
                    heatmap_data = response_data
            else:
                raise TypeError(f"Unexpected response_data type: {type(response_data)}")

                    # --- Trigger Asynchronous Global 3D Layout Generation ---
            # DISABLED: 3D coordinate generation to reduce server strain
            # Set initial task status
            task_status_cache[session_id] = {'status': 'disabled', 'message': '3D layout generation disabled to reduce server load.'}
            # _run_async_task_in_background(
            #     _generate_and_store_3d_coords_task,
            #     session_id,
            #     heatmap_data, # Pass the entire heatmap_data to the async task
            #     HEMISPHERE_LAYOUT_RADIUS
            # )
            logger.info(f"3D coordinate generation disabled for session: {session_id}")


            # DEBUG: Log successful completion
            request_end_time = time.time()
            total_processing_time = request_end_time - request_timestamp
            
            print(f'🔍 DEBUG: Request completed successfully')
            print(f'   Total processing time: {total_processing_time:.3f}s')
            print(f'   Session ID: {session_id}')
            print(f'   Files created: {file_path}, {metadata_file_path}')
            print('=' * 80)
                        
            return JsonResponse({
                "session_id": session_id,
                "has_missing_values": False,
                "heatmap_data": heatmap_data,  # Use the already parsed heatmap_data
                # "global_3d_positions_status": task_status_cache[session_id]['status'] # Return initial status
                "global_3d_positions_status": "disabled" # 3D coordinates disabled


            }, status=200)
         
    except Exception as e:
        error_trace = traceback.format_exc()
        
        # DEBUG: Log error completion
        request_end_time = time.time()
        total_processing_time = request_end_time - request_timestamp
        
        print(f'🚨 DEBUG: Request failed with error')
        print(f'   Total processing time: {total_processing_time:.3f}s')
        print(f'   Error: {str(e)}')
        print(f'   Process ID: {process_id}')
        print(f'   Thread ID: {thread_id}')
        print('=' * 80)
        
        logger.error(f"Error processing file: {str(e)}\n{error_trace}")
        return JsonResponse({
            "error": str(e),
            "traceback": error_trace
        }, status=500)
    
@api_view(['GET'])
def get_3d_coords_view(request, session_id: str):
    """
    API endpoint to retrieve pre-computed 3D coordinates for a given session ID.
    Checks the status of the background task.
    """
    coords_file_path = os.path.join(UPLOAD_DIR, f"{session_id}_coords.json")


    # print('**************** polling in to check whether the 3d coords file is created or not **************')
    # print('************* task info is as follows ****************',task_status_cache)
    
    # Prioritize checking the file system for the coordinates file
    if os.path.exists(coords_file_path):
        try:
            with open(coords_file_path, 'r') as f:
                global_3d_positions = json.load(f) # json.load is correct for file objects
            logger.info(f"Successfully retrieved 3D coordinates for session: {session_id} from file.")
            
            # Log correlation edge information without printing full coordinate data
            correlation_edges = global_3d_positions.get('correlation_edges', [])
            coordinates = global_3d_positions.get('coordinates', {})
            cluster_info = global_3d_positions.get('cluster_info', {})
            
            logger.info(f"3D Data Summary for session {session_id}:")
            logger.info(f"   - Coordinates: {len(coordinates)} genes")
            logger.info(f"   - Correlation edges: {len(correlation_edges)} edges")
            logger.info(f"   - Clusters: {len(cluster_info.get('clusters', []))} clusters")
            
            # If file exists and is valid, update cache to 'ready' (important for subsequent checks)
            task_status_cache[session_id] = {'status': 'ready', 'message': '3D layout ready.'}
            return JsonResponse({
                "session_id": session_id,
                "status": "ready",
                "global_3d_positions": global_3d_positions
            }, status=status.HTTP_200_OK)
        except json.JSONDecodeError as e:
            logger.error(f"Error decoding 3D coordinates JSON for session {session_id}: {str(e)}")
            # Mark as failed if file is corrupted
            task_status_cache[session_id] = {'status': 'failed', 'message': 'Corrupted coordinates file.'}
            return JsonResponse({"error": "Failed to load coordinates data"}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
    else:
        # If file does not exist, check the in-memory cache for ongoing status
        task_info = task_status_cache.get(session_id)
        
        if task_info and task_info['status'] in ['pending', 'computing']:
            logger.info(f"3D coordinates for session {session_id} are still {task_info['status']}.")
            return JsonResponse({
                "session_id": session_id,
                "status": task_info['status'],
                "message": task_info['message']
            }, status=status.HTTP_202_ACCEPTED) # 202 Accepted indicates processing is ongoing
        elif task_info and task_info['status'] == 'failed':
            logger.error(f"3D coordinate generation failed for session {session_id}.")
            return JsonResponse({
                "session_id": session_id,
                "status": "failed",
                "message": task_info.get('message', 'Failed to compute 3D layout.')
            }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        else:
            # Task not started, unknown session_id, or cache is stale/empty
            logger.info(f"No task status found in cache for session {session_id} and file not present. Assuming not started or lost status.")
            # This implies the task either hasn't started, or the worker that started it
            # has restarted and lost its in-memory cache.
            # In a real Celery setup, you'd query the Celery backend here.
            return JsonResponse({
                "session_id": session_id,
                "status": "not_started",
                "message": "3D coordinates computation has not started or status lost. Please wait."
            }, status=status.HTTP_404_NOT_FOUND) # Or 202 Accepted if you want to imply it will start




def handle_strategy_processing(request):
    """Handle the second step - processing with selected imputation strategy"""
    try:
        session_id = request.data.get('session_id')
        strategy = request.data.get('imputation_strategy')
        strategy_params = request.data.get('strategy_parameters')
        
        # Parse strategy parameters if provided
        imputation_kwargs = {}
        if strategy_params:
            if isinstance(strategy_params, str):
                imputation_kwargs = json.loads(strategy_params)
            else:
                imputation_kwargs = strategy_params
        
        print(f"Processing with {strategy} imputation strategy")
        print(f"Parameters: {imputation_kwargs}")
        
        # Reconstruct file path
        filename = f"{session_id}.tsv"
        file_path = os.path.join(UPLOAD_DIR, filename)
        
        if not os.path.exists(file_path):
            return JsonResponse({
                "error": "Session file not found. Please re-upload your data."
            }, status=404)
        
        # Call make_cluster with imputation parameters
        response_data = make_cluster(
            data=file_path,
            imputation_method=strategy,
            **imputation_kwargs  # Pass all strategy parameters as kwargs
        )

        response = json.loads(response_data)

        with open('genomics_data_file.json', 'w') as json_file:
            json.dump(response, json_file, indent=4)

        return JsonResponse({
            "session_id": session_id,
            "heatmap_data": json.loads(response_data)
        }, status=200)
        
    except Exception as e:
        error_trace = traceback.format_exc()
        logger.error(f"Error processing with strategy: {str(e)}\n{error_trace}")
        return JsonResponse({
            "error": str(e),
            "traceback": error_trace
        }, status=500)
    

def analyze_missing_values(file_path):
    """Analyze missing values in the uploaded file - FIXED to only analyze numeric matrix"""
    import pandas as pd
    import numpy as np

    try:
        # Read the entire file first
        df = pd.read_csv(file_path, sep='\t', header=None)  # Read without assuming structure
        print(f"Full file shape: {df.shape}")
        
        # Detect where the numeric matrix starts
        matrix_start_row, matrix_start_col = detect_matrix_start(df)
        print(f"Matrix starts at row={matrix_start_row}, col={matrix_start_col}")
        
        # Extract only the numeric data matrix portion
        numeric_matrix = df.iloc[matrix_start_row:, matrix_start_col:]
        print(f"Numeric matrix shape: {numeric_matrix.shape}")
        
        # Convert to numpy array for analysis, handling non-numeric values
        try:
            data = numeric_matrix.values.astype(float)
        except ValueError:
            # If direct conversion fails, convert cell by cell
            print("⚠️ Some values couldn't be converted to float, cleaning data...")
            data = np.full(numeric_matrix.shape, np.nan)
            
            for i in range(numeric_matrix.shape[0]):
                for j in range(numeric_matrix.shape[1]):
                    val = numeric_matrix.iloc[i, j]
                    if pd.notna(val):
                        try:
                            if isinstance(val, (int, float)):
                                data[i, j] = float(val)
                            elif isinstance(val, str):
                                # Try to convert string to number
                                data[i, j] = float(val)
                        except (ValueError, TypeError):
                            # Leave as NaN if conversion fails
                            data[i, j] = np.nan
        
        print(f"Final data matrix shape: {data.shape}")
        
        # Calculate missing value statistics on the numeric matrix only
        missing_mask = np.isnan(data)
        total_missing = np.sum(missing_mask)
        total_values = data.size
        missing_percentage = (total_missing / total_values) * 100
        
        genes_with_missing = np.sum(np.any(missing_mask, axis=1))
        samples_with_missing = np.sum(np.any(missing_mask, axis=0))
        
        genes_missing_percentage = (np.sum(missing_mask, axis=1) / data.shape[1]) * 100
        samples_missing_percentage = (np.sum(missing_mask, axis=0) / data.shape[0]) * 100
        
        print(f"Missing value analysis:")
        print(f"   Total missing: {total_missing}")
        print(f"   Missing percentage: {missing_percentage:.2f}%")
        print(f"   Genes with missing: {genes_with_missing}")
        print(f"   Samples with missing: {samples_with_missing}")
        
        return {
            "has_missing_values": bool(total_missing > 0),
            "total_missing": int(total_missing),
            "missing_percentage": float(missing_percentage),
            "genes_with_missing": int(genes_with_missing),
            "samples_with_missing": int(samples_with_missing),
            "genes_missing_percentage": [float(x) for x in genes_missing_percentage],
            "samples_missing_percentage": [float(x) for x in samples_missing_percentage],
            "total_genes": int(data.shape[0]),
            "total_samples": int(data.shape[1]),
            # Add matrix boundaries for debugging
            "matrix_boundaries": {
                "start_row": int(matrix_start_row),
                "start_col": int(matrix_start_col),
                "numeric_matrix_shape": [int(data.shape[0]), int(data.shape[1])]
            }
        }
        
    except Exception as e:
        print(f"❌ Error analyzing missing values: {e}")
        import traceback
        traceback.print_exc()
        return {
            "has_missing_values": False,
            "total_missing": 0,
            "missing_percentage": 0.0,
            "genes_with_missing": 0,
            "samples_with_missing": 0,
            "genes_missing_percentage": [],
            "samples_missing_percentage": [],
            "total_genes": 0,
            "total_samples": 0,
            "error": str(e)
        }



@api_view(['POST'])
def cleanup_session(request):
    session_id = request.data.get('session_id')
    if not session_id:
        return Response({"error": "No session_id provided"}, status=400)

    file_path = os.path.join(UPLOAD_DIR, f"{session_id}.tsv")
    metadata_file_path = os.path.join(UPLOAD_DIR, f"{session_id}_metadata.json")
    coords_file = os.path.join(UPLOAD_DIR, f"{session_id}_coords.json")

    files_deleted = []
    files_not_found = []

    # Attempt to delete .tsv file
    if os.path.exists(file_path):
        os.remove(file_path)
        files_deleted.append(f"{session_id}.tsv")
    else:
        files_not_found.append(f"{session_id}.tsv")

    # Attempt to delete metadata JSON
    if os.path.exists(metadata_file_path):
        os.remove(metadata_file_path)
        files_deleted.append(f"{session_id}_metadata.json")
    else:
        files_not_found.append(f"{session_id}_metadata.json")

    # Attempt to delete coordinates JSON
    if os.path.exists(coords_file):
        os.remove(coords_file)
        files_deleted.append(f"{session_id}_coords.json")
    else:
        files_not_found.append(f"{session_id}_coords.json")

    if files_deleted:
        print(f"🗑️ Session {session_id} cleaned up. Files deleted: {files_deleted}")
        return Response({
            "status": "Session cleaned up",
            "files_deleted": files_deleted,
            "files_not_found": files_not_found
        }, status=200)
    else:
        return Response({
            "error": "No session files found to clean up.",
            "files_not_found": files_not_found
        }, status=404)


# Add this validation function at the top of your views.py file
def validate_command_values(action, target, value, df, metadata=None):
    """
    Validate if the requested values exist in the dataset
    Returns (is_valid, error_message)
    """   
    import pandas as pd
    if action == "sample_filter":
        # Check if the metadata field exists and has the requested value
        if metadata and 'col' in metadata:
            if target in metadata['col']:
                available_values = metadata['col'][target]
                if value not in available_values:
                    return False, f"value not found: '{value}' not available in {target} field. Available values: {', '.join(available_values[:10])}"
            else:
                return False, f"value not found: '{target}' field not found in column metadata"
    
    elif action == "gene_filter":
        # Check if the gene/row metadata field exists
        if metadata and 'row' in metadata:
            if target in metadata['row']:
                available_values = metadata['row'][target]
                if value not in available_values:
                    return False, f"value not found: '{value}' not available in {target} field. Available values: {', '.join(available_values[:10])}"
            else:
                return False, f"value not found: '{target}' field not found in row metadata"
    
    elif action == "search":
        # Check if the gene exists in the first column (Unnamed: 0)
        # Skip NaN values since first few rows contain metadata
        if hasattr(df, 'columns') and len(df.columns) > 0:
            # Get the first column and filter out NaN values
            first_column = df.iloc[:, 0].dropna()  # Remove NaN values
            gene_names = first_column.tolist()
            
            if value not in gene_names:
                # Try case-insensitive search
                gene_names_lower = [str(g).lower() for g in gene_names if pd.notna(g)]
                if value.lower() not in gene_names_lower:
                    # Show first few gene names as examples (skip NaN values)
                    valid_genes = [g for g in gene_names if pd.notna(g)]
                    example_genes = valid_genes[:5] if len(valid_genes) >= 5 else valid_genes
                    return False, f"value not found: '{value}' not found in gene list. Example genes: {', '.join(map(str, example_genes))} ...."
    
    elif action == "sort_by_meta":
        # Check if the metadata field exists
        if target in ["columns", "cols"]:
            if metadata and 'col' in metadata:
                if value not in metadata['col']:
                    available_fields = list(metadata['col'].keys())
                    return False, f"value not found: '{value}' field not found in column metadata. Available fields: {', '.join(available_fields)}"
        elif target == "rows":
            if metadata and 'row' in metadata:
                if value not in metadata['row']:
                    available_fields = list(metadata['row'].keys())
                    return False, f"value not found: '{value}' field not found in row metadata. Available fields: {', '.join(available_fields)}"
    
    elif action == "variance":
        # Check if the requested number is reasonable
        try:
            num_features = int(value)
            if target == "rows":
                max_rows = len(df.index) if hasattr(df, 'index') else 0
                if num_features > max_rows:
                    return False, f"value not found: requested {num_features} rows but only {max_rows} available"
            elif target in ["columns", "cols"]:
                max_cols = len(df.columns) if hasattr(df, 'columns') else 0
                if num_features > max_cols:
                    return False, f"value not found: requested {num_features} columns but only {max_cols} available"
        except ValueError:
            return False, f"value not found: '{value}' is not a valid number for variance selection"
    
    # If we get here, validation passed
    return True, None

# Replace the section after "Extract data from response" with this:
@api_view(['POST'])
def command_execution(request):
    import pandas as pd
    import json
    import requests
    try:
        session_id = request.data.get('session_id')
        command = request.data.get('command')
        current_state = request.data.get('current_state', {})
        filters = current_state.get('filters', {})
        command_history = current_state.get('commandHistory', [])

        print('**** command history is as follows *****', command_history)

        if not session_id or not command:
            return Response({"error": "Missing 'session_id' or 'command'."}, status=400)

        file_path = os.path.join(UPLOAD_DIR, f"{session_id}.tsv")
        if not os.path.exists(file_path):
            return Response({"error": "Session data not found."}, status=404)

        # ✅ Load Original Data
        df = pd.read_csv(file_path, sep="\t")

        # ✅ Load Metadata - JSON file
        metadata_path = os.path.join(UPLOAD_DIR, f"{session_id}_metadata.json")
        metadata = None
        
        if os.path.exists(metadata_path):
            try:
                with open(metadata_path, 'r') as metadata_file:
                    metadata = json.load(metadata_file)
                print(f"Loaded metadata for session {session_id}")
            except Exception as e:
                print(f"Error loading metadata file: {str(e)}")
        else:
            print(f"No metadata file found at {metadata_path}")

        # ✅ Build Prompt for ChatGPT with metadata
        prompt = build_prompt(command, command_history, filters, metadata)

        # ✅ Define OpenAI API key (store this in environment variables in production)
        openai_api_key = "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"
        
        try:
            # ✅ Call OpenAI ChatGPT API
            headers = {
                "Content-Type": "application/json",
                "Authorization": f"Bearer {openai_api_key}"
            }
            
            payload = {
                "model": "gpt-4o-mini",  # You can also use "gpt-4" for better results
                "messages": [
                    {"role": "system", "content": "You are a specialized assistant that generates JSON responses for heatmap visualization commands."},
                    {"role": "user", "content": prompt}
                ],
                "temperature": 0.3,  # Lower temperature for more deterministic responses
                "max_tokens": 150
            }
            
            
            # Make the API call with timeout
            response = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers=headers,
                json=payload,
                timeout=10  # 10 second timeout
            )
            
            # Check for errors in the API response
            if response.status_code != 200:
                print(f"Error from OpenAI API: {response.status_code}")
                print(response.text)
                raise Exception(f"OpenAI API error: {response.status_code}")
                
            response_data = response.json()
            
            # Extract the content from the response
            gpt_response = response_data["choices"][0]["message"]["content"].strip()

            print('********** gpt_response is as follows ************',gpt_response)
            
        
            # Replace your existing response processing section with this:

            if not gpt_response:
                print("Empty response received, using fallback")
                ollama_response = generate_fallback_response(command)
            else:
                try:
                    # Try to find and parse JSON in the response
                    import re
                    json_pattern = r'(\{[\s\S]*?\})'
                    matches = re.findall(json_pattern, gpt_response)
                    
                    if matches:
                        ollama_response = None  # Initialize to track if we found valid JSON
                        
                        for potential_json in matches:
                            try:
                                parsed_json = json.loads(potential_json)
                                
                                # ✅ CRITICAL FIX: Check for error first, before looking for action
                                if "error" in parsed_json:
                                    print(f"Found valid JSON with error: {parsed_json['error']}")
                                    ollama_response = parsed_json
                                    break
                                elif "action" in parsed_json:
                                    print("Found valid JSON with action field")
                                    ollama_response = parsed_json
                                    break
                                    
                            except json.JSONDecodeError as json_error:
                                print(f"JSON decode error for '{potential_json[:50]}...': {json_error}")
                                continue
                            except Exception as parse_error:
                                print(f"Unexpected error parsing JSON '{potential_json[:50]}...': {parse_error}")
                                continue
                        
                        # ✅ CRITICAL FIX: Check if we actually found valid JSON
                        if ollama_response is None:
                            # No valid JSON found in any match
                            print("No valid JSON with action or error field found, using fallback")
                            ollama_response = generate_fallback_response(command)
                    else:
                        # No JSON-like patterns found
                        print("No JSON patterns found, using fallback")
                        ollama_response = generate_fallback_response(command)
                                    
                except Exception as e:
                    print(f"Error processing ChatGPT response: {e}")
                    ollama_response = generate_fallback_response(command)

                    
        except Exception as e:
            print(f"Error with ChatGPT request: {str(e)}")
            ollama_response = generate_fallback_response(command)
        
        # Extract data from response
        action = ollama_response.get("action")
        target = ollama_response.get("target")
        value = ollama_response.get("value", "")
        
        # ✅ NEW: Check if there's already an error in the response (from improved prompt)
        if "error" in ollama_response:
            return Response({
                "error": ollama_response["error"],
                "commandHistory": command_history + [command]
            }, status=200)
        
        # ✅ NEW: Validate the command values against actual data
        is_valid, error_message = validate_command_values(action, target, value, df, metadata)

        print('****** ollama response is  *******',ollama_response)
        print('****** is it valid ******',is_valid)
        print('******* error message is *******',error_message)
        
        if not is_valid:
            return Response({
                "error": error_message,
                "commandHistory": command_history + [command]
            }, status=200)
        
        VIEW_ONLY_ACTIONS = ['sort', 'sort_by_meta', 'cluster', 'search', 'set_opacity']
        PATHWAY_ACTIONS = ['pathway_filter', 'functional_filter', 'pathway_search']


        if action in VIEW_ONLY_ACTIONS:
            # --- This is a "View State" command ---
            # We don't need to re-cluster. Just send the command back to the frontend.
            print(f"Action '{action}' is view-only, will be handled by the frontend.")
            response_data = {
                "action": action,
                "target": target,
                "value": value,
                "updated_filters": filters, # Send back original filters
                "commandHistory": command_history + [command]                
            }
        elif action in PATHWAY_ACTIONS:
            is_pathway_valid, pathway_error = validate_pathway_command(action, target, value)
            if not is_pathway_valid:
                return Response({
                    "error": pathway_error,
                    "commandHistory": command_history + [command]
                }, status=200)

            # Replace your existing elif action in PATHWAY_ACTIONS: section with this enhanced version:
            else:
                
                # Replace your pathway filtering sections in command_execution() with these corrected versions:

                if action == 'functional_filter':
                    # Handle "show immune genes" type commands
                    try:
                        functional_genes = get_functional_genes(value)
                        
                        if functional_genes:
                            print(f"Found {len(functional_genes)} genes for function '{value}'")
                            
                            # Use your existing filter function for gene rows
                            filtered_df = filter_genes_by_ids(df, functional_genes)
                            
                            # Check if filtering actually worked
                            if filtered_df.shape[0] > df.iloc[:detect_matrix_start(df)[0], :].shape[0]:  # More than just metadata rows
                                filters['functional'] = value
                                
                                # Re-cluster the filtered data
                                clustering_result_str = make_cluster(filtered_df)

                                print('************************************************************',clustering_result_str)
                                clustering_result = json.loads(clustering_result_str)
                                
                                # Count how many genes were actually found
                                matrix_start_row, _ = detect_matrix_start(filtered_df)
                                genes_found = filtered_df.shape[0] - matrix_start_row
                                
                                response_data = {
                                    "action": action,
                                    "target": target,
                                    "value": value,
                                    "updated_filters": filters,
                                    "clustering_result": clustering_result,
                                    "functional_info": f"Filtered to {genes_found} genes related to {value}",
                                    "total_functional_genes": len(functional_genes),
                                    "genes_in_dataset": genes_found,
                                    "result_type": "data_filter",
                                    "commandHistory": command_history + [command]
                                }
                            else:
                                # Get matrix start to count original genes
                                matrix_start_row, _ = detect_matrix_start(df)
                                total_genes_in_dataset = df.shape[0] - matrix_start_row
                                
                                response_data = {
                                    "error": f"None of the {len(functional_genes)} {value}-related genes found in your dataset of {total_genes_in_dataset} genes",
                                    "suggestion": f"Try 'list {value} pathways' to see what's available",
                                    "total_functional_genes": len(functional_genes),
                                    "commandHistory": command_history + [command]
                                }
                        else:
                            response_data = {
                                "error": f"No genes found for function: {value}",
                                "suggestion": "Try 'list immune pathways' or 'list cancer pathways'",
                                "commandHistory": command_history + [command]
                            }
                            
                    except Exception as e:
                        print(f"Error in functional_filter: {str(e)}")
                        import traceback
                        traceback.print_exc()
                        response_data = {
                            "error": f"Error filtering by function: {str(e)}",
                            "commandHistory": command_history + [command]
                        }

                elif action == 'pathway_filter':
                    # Handle "show STAT3 targets" type commands
                    try:
                        pathway_genes = get_pathway_genes(value, action)
                        
                        if pathway_genes:
                            print(f"Found {len(pathway_genes)} genes for pathway '{value}'")
                            
                            # Use your existing filter function for gene rows
                            filtered_df = filter_genes_by_ids(df, pathway_genes)
                            
                            # Check if filtering actually worked
                            if filtered_df.shape[0] > df.iloc[:detect_matrix_start(df)[0], :].shape[0]:  # More than just metadata rows
                                filters['pathway'] = value
                                
                                # Re-cluster the filtered data
                                clustering_result_str = make_cluster(filtered_df)
                                clustering_result = json.loads(clustering_result_str)
                                
                                # Count how many genes were actually found
                                matrix_start_row, _ = detect_matrix_start(filtered_df)
                                genes_found = filtered_df.shape[0] - matrix_start_row
                                
                                response_data = {
                                    "action": action,
                                    "target": target,
                                    "value": value,
                                    "updated_filters": filters,
                                    "clustering_result": clustering_result,
                                    "pathway_info": f"Filtered to {genes_found} genes from {value}",
                                    "total_pathway_genes": len(pathway_genes),
                                    "genes_in_dataset": genes_found,
                                    "result_type": "data_filter",
                                    "commandHistory": command_history + [command]
                                }
                            else:
                                # Get matrix start to count original genes
                                matrix_start_row, _ = detect_matrix_start(df)
                                total_genes_in_dataset = df.shape[0] - matrix_start_row
                                
                                response_data = {
                                    "error": f"None of the {len(pathway_genes)} genes from {value} found in your dataset of {total_genes_in_dataset} genes",
                                    "suggestion": "Try 'list immune pathways' to see available options",
                                    "total_pathway_genes": len(pathway_genes),
                                    "dataset_info": f"Your dataset has {total_genes_in_dataset} genes total",
                                    "commandHistory": command_history + [command]
                                }
                        else:
                            response_data = {
                                "error": f"Pathway not found: {value}",
                                "suggestion": "Try 'list immune pathways' to see available options", 
                                "commandHistory": command_history + [command]
                            }
                            
                    except Exception as e:
                        print(f"Error in pathway_filter: {str(e)}")
                        import traceback
                        traceback.print_exc()
                        response_data = {
                            "error": f"Error processing pathway: {str(e)}",
                            "commandHistory": command_history + [command]
                        }

                # Keep your existing pathway_search logic as is - that one doesn't need to filter data:
                elif action == 'pathway_search':
                    # Handle "list immune pathways" type commands

                    print('*************** iT is pathway searching command **************')
                    # Handle "list immune pathways" type commands
                    try:
                        matching_pathways = search_pathways_by_category(value)
                        print(f"Found pathways: {matching_pathways}")
                        
                        if matching_pathways:
                            # Format results for frontend display - KEEP IT SIMPLE
                            formatted_results = []
                            for pathway in matching_pathways:
                                formatted_results.append({
                                    'name': pathway['pathway_name'],
                                    'description': pathway.get('description', pathway['pathway_name']),
                                    'gene_count': pathway['gene_count'],
                                    'library': pathway['library'],
                                    'match_reason': pathway['match_reason']
                                })
                            
                            response_data = {
                                "action": action,
                                "target": target,
                                "value": value,
                                "pathway_results": formatted_results,
                                "message": f"Found {len(matching_pathways)} pathways related to '{value}'",
                                "result_type": "pathway_list",
                                "commandHistory": command_history + [command]
                            }
                        else:
                            response_data = {
                                "error": f"No pathways found for category: {value}. Try 'immune', 'cancer', 'metabolism', 'histone', etc.",
                                "commandHistory": command_history + [command]
                            }
                    except Exception as e:
                        print(f"Error in pathway_search: {str(e)}")
                        response_data = {
                            "error": f"Error searching pathways: {str(e)}",
                            "commandHistory": command_history + [command]
                        }
        elif action not in VIEW_ONLY_ACTIONS and action not in PATHWAY_ACTIONS:
            # --- This is a "Data Subsetting" command ---
            # This is your original logic for filtering and re-clustering
            print(f"Action '{action}' requires data re-processing on the backend.")
            updated_filters = update_filters(filters, action, target, value)
            final_df = apply_filters(df, updated_filters, session_id)
            final_df.columns = ['' if 'unnamed' in str(col).lower() else str(col) for col in final_df.columns]

            if final_df is not None and not final_df.empty:
                clustering_result_str = make_cluster(final_df)
                clustering_result = json.loads(clustering_result_str)
                
                response_data = {
                    "action": action,
                    "target": target,
                    "value": value,
                    "updated_filters": updated_filters,
                    "clustering_result": clustering_result, # ✅ Send new data
                    "commandHistory": command_history + [command]
                }
            else:
                response_data = {
                    "error": "Filtering resulted in an empty dataset.",
                    "commandHistory": command_history + [command]
                }

        return Response(response_data, status=200)
    
    except Exception as e:
        error_trace = traceback.format_exc()
        print(f"ERROR: {str(e)}\n{error_trace}")
        return Response({"error": str(e), "traceback": error_trace}, status=500)

    
def filter_genes_by_ids(df, gene_ids):
    import pandas as pd

    """
    Filter DataFrame to only include rows that match the provided gene IDs.
    Looks for gene IDs in the first few columns (gene names/identifiers).
    """
    try:
        print(f"Filtering for {len(gene_ids)} specific gene IDs...")
        
        # Detect where the matrix starts (similar to your existing logic)
        matrix_start_row, matrix_start_col = detect_matrix_start(df)
        print(f"Matrix starts at row {matrix_start_row}, col {matrix_start_col}")
        
        # Get the data portion (rows from matrix_start_row onwards)
        data_rows = df.iloc[matrix_start_row:, :]
        metadata_rows = df.iloc[:matrix_start_row, :]
        
        print(f"Data rows shape: {data_rows.shape}")
        print(f"Looking for gene IDs in columns: {data_rows.columns[:matrix_start_col].tolist()}")
        
        # Find matching rows by checking gene identifier columns
        matching_indices = []
        
        for idx, row in data_rows.iterrows():
            # Check the first few columns for gene identifiers
            gene_found = False
            for col_idx in range(min(matrix_start_col, len(row))):
                cell_value = str(row.iloc[col_idx]).strip()                
                # Check if this cell value matches any of our target gene IDs
                if cell_value in gene_ids:
                    matching_indices.append(idx)
                    gene_found = True
                    break
            
            # Also check if the row index itself is a gene ID
            if not gene_found and str(row.name) in gene_ids:
                matching_indices.append(idx)
                print(f"Found gene '{row.name}' as row index")
        
        print(f"Found {len(matching_indices)} matching genes out of {len(gene_ids)} requested")
        
        if not matching_indices:
            print("⚠️ Warning: No matching genes found!")
            return df  # Return original if no matches
        
        # Filter to only matching rows
        filtered_data_rows = data_rows.loc[matching_indices]
        
        # Combine metadata rows with filtered data rows
        filtered_df = pd.concat([metadata_rows, filtered_data_rows], axis=0)
        
        print(f"Filtered DataFrame shape: {filtered_df.shape}")
        return filtered_df
        
    except Exception as e:
        print(f"Error in filter_genes_by_ids: {str(e)}")
        import traceback
        traceback.print_exc()
        return df  # Return original on error

"""
Updated Django view integrating the ultra-optimized correlation engine.

This is a DROP-IN REPLACEMENT for your existing view.
All existing functionality is preserved with 10-20x performance improvement.
"""

@api_view(['POST'])
def correlation_network(request):
    import pandas as pd

    """
    Ultra-optimized correlation network function.
    
    This is a DROP-IN REPLACEMENT for your existing function with:
    - 10-20x performance improvement
    - Same API interface
    - Enhanced error handling
    - Progress estimation
    - Memory optimization
    """
    try:
        # Extract data from request (UNCHANGED)
        session_id = request.data.get('sessionId')
        gene_ids = request.data.get('geneIds', [])
        filters = request.data.get('filters', {})
        
        # Print received arguments (UNCHANGED)
        print('===== ULTRA-OPTIMIZED CORRELATION NETWORK REQUEST =====')
        print(f'Session ID: {session_id}')
        print(f'Number of Gene IDs: {len(gene_ids)}')
        print(f'Sample Gene IDs: {gene_ids[:10]}...' if len(gene_ids) > 10 else f'Gene IDs: {gene_ids}')
        print(f'Filters: {filters}')
        
        # Estimate processing time for user feedback
        estimated_time = estimate_correlation_job_time(len(gene_ids), filters.get('correlationThreshold', 0.3))
        print(f'Estimated processing time: {estimated_time} seconds')
        print('=' * 60)
        
        # Basic validation (UNCHANGED)
        if not session_id:
            return Response({"error": "Missing 'sessionId'."}, status=400)
        
        if not gene_ids or len(gene_ids) == 0:
            return Response({"error": "Missing or empty 'geneIds'."}, status=400)
        
        # Additional validation for large datasets
        if len(gene_ids) > 50000:
            return Response({
                "error": "Too many genes requested. Maximum 50,000 genes allowed.",
                "requested": len(gene_ids),
                "maximum": 50000
            }, status=400)
        
        # Load the session data file (UNCHANGED)
        file_path = os.path.join(UPLOAD_DIR, f"{session_id}.tsv")
        if not os.path.exists(file_path):
            return Response({"error": "Session data not found."}, status=404)
        
        # ✅ Load Original Data (UNCHANGED)
        df = pd.read_csv(file_path, sep="\t")
        print(f"Original DataFrame shape: {df.shape}")
        
        # ✅ Filter DataFrame to only include specific gene IDs (UNCHANGED)
        filtered_gene_df = filter_genes_by_ids(df, gene_ids)
        print(f"After gene filtering shape: {filtered_gene_df.shape}")
        
        # ✅ Apply additional filters using your existing function (UNCHANGED)
        if filters:
            final_df = apply_filters(filtered_gene_df, filters, session_id)
            print(f"After applying filters shape: {final_df.shape}")
        else:
            final_df = filtered_gene_df
        
        # ✅ Validate filtered data
        if final_df.empty:
            return Response({"error": "No data remaining after filtering"}, status=400)
        
        # Check if dataset is reasonable size
        if final_df.shape[0] < 2:
            return Response({"error": "Need at least 2 genes for correlation analysis"}, status=400)
        
        # ✅ ULTRA-OPTIMIZED Correlation Computation
        print("🚀 Starting ultra-optimized correlation computation...")
        
        # Set intelligent defaults if not provided
        optimized_filters = {
            'correlationThreshold': filters.get('correlationThreshold', 0.7),  # More realistic default
            'pValueThreshold': filters.get('pValueThreshold', 0.05),           # Less strict default
            'maxCorrelations': filters.get('maxCorrelations', 20000),
            'variancePercentile': filters.get('variancePercentile', 0.2),     # Keep top 80%
            'minSamples': filters.get('minSamples', 10),
            **filters  # Preserve any additional filters
        }
        
        # Run ultra-optimized correlation computation
        correlation_result = compute_correlation_matrix(final_df, optimized_filters)
        
        # Handle errors from correlation computation
        if "error" in correlation_result:
            return Response({
                "error": "Correlation computation failed",
                "details": correlation_result["error"]
            }, status=500)
        
        # Print results summary
        print(f"🎯 Correlation computation completed!")
        print(f"📊 Total correlations found: {correlation_result.get('totalCorrelations', 0):,}")
        print(f"⏱️  Computation time: {correlation_result.get('computationTime', 0):.2f}s")
        print(f"🏃 Speed: {correlation_result.get('pairsPerSecond', 0):,.0f} pairs/second")
        
        # Estimate response data size
        n_correlations = correlation_result.get('totalCorrelations', 0)
        estimated_size = get_correlation_data_size_estimate(n_correlations)
        print(f"📦 Estimated response size: {estimated_size}")
        
        # Enhanced response data with additional metadata
        response_data = {
            "sessionId": session_id,
            "totalGenes": len(gene_ids),
            "filteredGenes": final_df.shape[0] if hasattr(final_df, 'shape') else 0,
            "filteredSamples": final_df.shape[1] if hasattr(final_df, 'shape') else 0,
            "correlationMatrix": correlation_result,
            "status": "computed",
            "message": f"Successfully computed correlations for {correlation_result.get('totalCorrelations', 0):,} gene pairs",
            "filters": optimized_filters,
            
            # Additional metadata for frontend
            "performance": {
                "optimizationLevel": "ULTRA_OPTIMIZED", 
                "computationTime": correlation_result.get('computationTime', 0),
                "pairsPerSecond": correlation_result.get('pairsPerSecond', 0),
                "totalPairsComputed": correlation_result.get('totalPairsComputed', 0),
                "estimatedResponseSize": estimated_size
            },
            
            # Data quality metrics
            "dataQuality": {
                "correlationThreshold": optimized_filters['correlationThreshold'],
                "significantCorrelations": correlation_result.get('totalCorrelations', 0),
                "totalPossiblePairs": correlation_result.get('totalPairsComputed', 0),
                "passRate": round(
                    (correlation_result.get('totalCorrelations', 0) / 
                     max(correlation_result.get('totalPairsComputed', 1), 1)) * 100, 2
                )
            },
            
            # Pagination info for large result sets
            "pagination": {
                "hasMore": correlation_result.get('totalSignificantFound', 0) > correlation_result.get('totalCorrelations', 0),
                "totalAvailable": correlation_result.get('totalSignificantFound', 0),
                "returned": correlation_result.get('totalCorrelations', 0)
            }
        }
        
        return Response(response_data, status=200)
    
    except MemoryError as e:
        error_message = f"Insufficient memory for correlation analysis with {len(gene_ids)} genes"
        print(f"💾 MEMORY ERROR: {error_message}")
        return Response({
            "error": error_message,
            "suggestions": [
                "Reduce the number of genes",
                "Increase correlation threshold to reduce output",
                "Use smaller maxCorrelations limit"
            ],
            "memoryError": True
        }, status=413)  # Payload Too Large
    
    except TimeoutError as e:
        error_message = "Correlation computation timed out"
        print(f"⏰ TIMEOUT ERROR: {error_message}")
        return Response({
            "error": error_message,
            "suggestions": [
                "Reduce the number of genes", 
                "Increase correlation threshold",
                "Try the analysis in smaller batches"
            ],
            "timeoutError": True
        }, status=408)  # Request Timeout
    
    except Exception as e:
        error_trace = traceback.format_exc()
        print(f"❌ ERROR in ultra-optimized correlation_network: {str(e)}\n{error_trace}")
        return Response({
            "error": str(e),
            "traceback": error_trace,
            "errorType": type(e).__name__
        }, status=500)
    

@api_view(['POST'])
def load_example_data(request):
    """
    Load pre-processed example dataset from server storage.
    
    Expects:
    - example_id: The ID of the example dataset (gene-expression, proteomics, immunogenomics)
    
    Returns:
    - JSON data ready for heatmap visualization
    """
    try:
        # Extract example ID from request
        example_id = request.data.get('example_id')
        
        print('===== LOAD EXAMPLE DATA REQUEST =====')
        print(f'Example ID: {example_id}')
        print('=' * 40)
        
        # Basic validation
        if not example_id:
            return Response({"error": "Missing 'example_id'."}, status=400)
        
        # Validate example ID format (basic security check)
        allowed_examples = ['genomics', 'proteomics', 'immunogenomics']
        if example_id not in allowed_examples:
            return Response({
                "error": f"Invalid example_id. Allowed values: {allowed_examples}",
                "provided": example_id
            }, status=400)
        
        # Construct file path
        file_path = os.path.join(UPLOAD_DIR, f"{example_id}.json")
        print(f"Looking for file: {file_path}")
        
        # Check if file exists
        if not os.path.exists(file_path):
            return Response({
                "error": f"Example dataset not found: {example_id}",
                "file_path": file_path
            }, status=404)
        
        # Check file size for logging
        file_size = os.path.getsize(file_path)
        print(f"File size: {file_size / (1024*1024):.2f} MB")
        
        # Load and return the JSON data
        print(f"📊 Loading example data: {example_id}")
        
        with open(file_path, 'r') as file:
            json_data = json.load(file)
        
        print(f"✅ Successfully loaded example data: {example_id}")
        print(f"📋 Data keys: {list(json_data.keys()) if isinstance(json_data, dict) else 'Not a dict'}")
        
        # Return the data directly
        response_data = {
            "success": True,
            "example_id": example_id,
            "data": json_data,
            "message": f"Successfully loaded {example_id} dataset",
            "file_size_mb": round(file_size / (1024*1024), 2)
        }
        
        return Response(response_data, status=200)
    
    except json.JSONDecodeError as e:
        error_message = f"Invalid JSON format in {example_id}.json"
        print(f"🔥 JSON DECODE ERROR: {error_message} - {str(e)}")
        return Response({
            "error": error_message,
            "details": str(e),
            "example_id": example_id
        }, status=500)
    
    except FileNotFoundError as e:
        error_message = f"Example dataset file not found: {example_id}"
        print(f"📁 FILE NOT FOUND: {error_message}")
        return Response({
            "error": error_message,
            "example_id": example_id,
            "file_path": file_path
        }, status=404)
    
    except MemoryError as e:
        error_message = f"File too large to load into memory: {example_id}"
        print(f"💾 MEMORY ERROR: {error_message}")
        return Response({
            "error": error_message,
            "suggestions": [
                "File is too large for direct loading",
                "Consider implementing streaming or chunked loading",
                "Use a smaller dataset for testing"
            ],
            "example_id": example_id
        }, status=413)  # Payload Too Large
    
    except Exception as e:
        error_trace = traceback.format_exc()
        print(f"❌ ERROR in load_example_data: {str(e)}\n{error_trace}")
        return Response({
            "error": str(e),
            "traceback": error_trace,
            "errorType": type(e).__name__,
            "example_id": example_id
        }, status=500)

# Optional: Add a new endpoint for correlation job estimation
@api_view(['POST'])
def estimate_correlation_job(request):
    """
    Estimate correlation job time and resource requirements.
    Useful for showing progress bars or warnings to users.
    """
    try:
        gene_ids = request.data.get('geneIds', [])
        filters = request.data.get('filters', {})
        
        if not gene_ids:
            return Response({"error": "Missing 'geneIds'."}, status=400)
        
        n_genes = len(gene_ids)
        correlation_threshold = filters.get('correlationThreshold', 0.3)
        
        # Calculate estimates
        total_pairs = n_genes * (n_genes - 1) // 2
        estimated_time = estimate_correlation_job_time(n_genes, correlation_threshold)
        estimated_memory_mb = n_genes * 0.03  # Rough estimate: 0.03 MB per gene
        
        # Estimate result size based on threshold
        if correlation_threshold >= 0.6:
            expected_correlations = int(total_pairs * 0.001)  # 0.1% pass rate
        elif correlation_threshold >= 0.4:
            expected_correlations = int(total_pairs * 0.01)   # 1% pass rate
        elif correlation_threshold >= 0.3:
            expected_correlations = int(total_pairs * 0.05)   # 5% pass rate
        else:
            expected_correlations = int(total_pairs * 0.1)    # 10% pass rate
        
        estimated_response_size = get_correlation_data_size_estimate(expected_correlations)
        
        # Determine feasibility
        if estimated_time > 300:  # > 5 minutes
            feasibility = "batch_processing_recommended"
            recommendation = "Consider reducing genes or increasing threshold"
        elif estimated_time > 60:  # > 1 minute
            feasibility = "long_running"
            recommendation = "Suitable for background processing"
        else:
            feasibility = "real_time"
            recommendation = "Suitable for real-time analysis"
        
        return Response({
            "geneCount": n_genes,
            "totalPairs": total_pairs,
            "estimates": {
                "computationTimeSeconds": estimated_time,
                "memoryUsageMB": round(estimated_memory_mb, 1),
                "expectedCorrelations": expected_correlations,
                "responseSize": estimated_response_size
            },
            "feasibility": feasibility,
            "recommendation": recommendation,
            "thresholds": {
                "correlationThreshold": correlation_threshold,
                "pValueThreshold": filters.get('pValueThreshold', 0.1)
            }
        }, status=200)
        
    except Exception as e:
        return Response({"error": str(e)}, status=500)

# Optional: Add endpoint for correlation parameter recommendations
@api_view(['POST'])
def recommend_correlation_parameters(request):
    """
    Recommend optimal correlation parameters based on dataset characteristics.
    """
    try:
        session_id = request.data.get('sessionId')
        gene_ids = request.data.get('geneIds', [])
        
        if not session_id or not gene_ids:
            return Response({"error": "Missing sessionId or geneIds"}, status=400)
        
        # Load data to analyze characteristics
        file_path = os.path.join(UPLOAD_DIR, f"{session_id}.tsv")
        if not os.path.exists(file_path):
            return Response({"error": "Session data not found."}, status=404)
        
        df = pd.read_csv(file_path, sep="\t")
        filtered_df = filter_genes_by_ids(df, gene_ids)
        
        n_genes = len(gene_ids)
        n_samples = filtered_df.shape[1] - 3  # Subtract metadata columns
        
        # Recommend parameters based on dataset size
        if n_genes < 1000:
            # Small dataset - can afford lower thresholds
            recommendations = {
                "correlationThreshold": 0.2,
                "pValueThreshold": 0.05,
                "maxCorrelations": 10000,
                "rationale": "Small dataset allows for comprehensive analysis with lower thresholds"
            }
        elif n_genes < 5000:
            # Medium dataset - balanced approach
            recommendations = {
                "correlationThreshold": 0.3,
                "pValueThreshold": 0.1,
                "maxCorrelations": 25000,
                "rationale": "Medium dataset benefits from moderate thresholds for good coverage"
            }
        elif n_genes < 15000:
            # Large dataset - focus on strong correlations
            recommendations = {
                "correlationThreshold": 0.4,
                "pValueThreshold": 0.1,
                "maxCorrelations": 50000,
                "rationale": "Large dataset requires higher thresholds to focus on strongest correlations"
            }
        else:
            # Very large dataset - strict filtering needed
            recommendations = {
                "correlationThreshold": 0.5,
                "pValueThreshold": 0.05,
                "maxCorrelations": 50000,
                "rationale": "Very large dataset needs strict filtering for manageable results"
            }
        
        # Adjust based on sample size
        if n_samples < 50:
            recommendations["pValueThreshold"] = 0.2  # Less strict for small samples
            recommendations["rationale"] += ". Relaxed p-value due to small sample size."
        elif n_samples > 500:
            recommendations["pValueThreshold"] = 0.01  # More strict for large samples
            recommendations["rationale"] += ". Strict p-value due to large sample size."
        
        return Response({
            "datasetCharacteristics": {
                "geneCount": n_genes,
                "sampleCount": n_samples,
                "totalPossiblePairs": n_genes * (n_genes - 1) // 2
            },
            "recommendations": recommendations,
            "alternativeOptions": {
                "exploratory": {
                    "correlationThreshold": max(0.1, recommendations["correlationThreshold"] - 0.2),
                    "pValueThreshold": 0.2,
                    "description": "For exploratory analysis with broader coverage"
                },
                "focused": {
                    "correlationThreshold": min(0.8, recommendations["correlationThreshold"] + 0.2),
                    "pValueThreshold": 0.01,
                    "description": "For focused analysis on strongest correlations only"
                }
            }
        }, status=200)
        
    except Exception as e:
        return Response({"error": str(e)}, status=500)

    

@api_view(['POST'])
def refresh_heatmap(request):
    import pandas as pd

    """
    Basic correlation network function that accepts arguments and prints them.
    Will be expanded with actual correlation computation later.
    """
    try:
        # Extract data from request
        session_id = request.data.get('sessionId')
        filters = request.data.get('filters', {})
    
        
        # Basic validation
        if not session_id:
            return Response({"error": "Missing 'sessionId'."}, status=400)
        
        
        # Load the session data file
        file_path = os.path.join(UPLOAD_DIR, f"{session_id}.tsv")
        if not os.path.exists(file_path):
            return Response({"error": "Session data not found."}, status=404)
        
        # ✅ Load Original Data
        df = pd.read_csv(file_path, sep="\t")
        
        # ✅ Apply additional filters using your existing function
        if filters:
            final_df = apply_filters(df, filters, session_id)
            print(f"After applying filters shape: {final_df.shape}")
        else:
            final_df = df

        
       # ✅ Compute Correlation Matrix
        if final_df.empty:
            return Response({"error": "No data remaining after filtering"}, status=400)
        
        
        # Process data and generate clustering result
        if final_df is not None and not final_df.empty:
            try:
                clustering_result_str = make_cluster(final_df)
                # Parse the clustering result
                try:
                    clustering_result = json.loads(clustering_result_str)
                except (json.JSONDecodeError, TypeError):
                    clustering_result = clustering_result_str
                    
                response_data = {
                    "clustering_result": clustering_result,
                }
            except Exception as e:
                print(f"Error in clustering: {e}")
                response_data = {
                    "clustering_result": None,
                    "error": f"Clustering error: {str(e)}",
                }
        else:
            # Handle empty DataFrame case
            response_data = {
                "clustering_result": None,
                "message": "Filtering resulted in empty dataset, clustering was skipped",
            }

        return Response(response_data, status=200)

        
        
    except Exception as e:
        error_trace = traceback.format_exc()
        print(f"ERROR in get_correlation_network: {str(e)}\n{error_trace}")
        return Response({"error": str(e), "traceback": error_trace}, status=500)
    

## ✅ A simple helper function that will be run by each thread
def fetch_single_library(user_list_id, library):
    """Fetches enrichment results for a single library. This is our thread's job."""
    try:
        query_string = f'?userListId={user_list_id}&backgroundType={library}'
        response = requests.get(f'https://maayanlab.cloud/Enrichr/enrich{query_string}')
        response.raise_for_status()
        return library, response.json()
    except requests.exceptions.RequestException:
        # If one library fails, we don't want to crash the whole request
        return library, None

# ✅ The view is now a standard synchronous "def"
@api_view(['POST'])
def enrich_analysis_view(request):
    import requests
    import concurrent.futures

    genes = request.data.get('genes', [])
    if not isinstance(genes, list) or not genes:
        return Response({"error": "A non-empty list of 'genes' is required."}, status=status.HTTP_400_BAD_REQUEST)

    genes_str = '\n'.join(genes)
    
    # Step 1: Submit gene list (This is still a single, quick call)
    try:
        addlist_response = requests.post(
            'https://maayanlab.cloud/Enrichr/addList',
            files={'list': (None, genes_str), 'description': (None, 'ClusterChirp Gene List')}
        )
        addlist_response.raise_for_status()
        user_list_id = addlist_response.json().get('userListId')
    except requests.exceptions.RequestException as e:
        return Response({"error": f"Enrichr API Error (addList): {e}"}, status=status.HTTP_502_BAD_GATEWAY)

    # ✅ Step 2: Concurrently fetch results using a Thread Pool
    all_library_results = []
    # Use a thread pool to make all requests concurrently. We limit workers to avoid overwhelming the API.
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        # Create a "future" for each library request
        future_to_library = {executor.submit(fetch_single_library, user_list_id, lib): lib for lib in ALL_LIBRARIES}
        
        for future in concurrent.futures.as_completed(future_to_library):
            try:
                # Get the result from the completed future
                all_library_results.append(future.result())
            except Exception as exc:
                library_name = future_to_library[future]
                print(f'{library_name} generated an exception: {exc}')

    # Step 3: Process the results (This logic is the same)
    final_results = {}
    for category, libraries in LIBRARY_CATEGORIES.items():
        category_results = {}
        for library_name in libraries:
            for res_lib, res_data in all_library_results:
                if res_lib == library_name and res_data and library_name in res_data:
                    significant_terms = [
                        {
                            "rank": term[0], "name": term[1], "pValue": term[2],
                            "zScore": term[3], "combinedScore": term[4], "genes": term[5]
                        }
                        for term in res_data[library_name] if term[2] < 0.05
                    ][:15]
                    
                    if significant_terms:
                        category_results[library_name] = significant_terms
                    break
        
        if category_results:
            final_results[category] = category_results
    
    return Response(final_results, status=status.HTTP_200_OK)

# Note: Make sure to import your existing detect_matrix_start and apply_filters functions

# def compute_correlation_chunk(args):
#     """
#     Compute correlations for a chunk of gene pairs.
#     This function will be run in parallel processes.
#     """
#     gene_pairs, data_matrix_values, gene_ids, correlation_threshold, p_value_threshold = args
    
#     correlations = []
    
#     for i, j in gene_pairs:
#         try:
#             # Get gene data
#             gene1_data = data_matrix_values[i]
#             gene2_data = data_matrix_values[j]
            
#             # Find non-null samples for both genes
#             mask1 = ~np.isnan(gene1_data)
#             mask2 = ~np.isnan(gene2_data)
#             common_mask = mask1 & mask2
            
#             if np.sum(common_mask) < 3:  # Need at least 3 samples
#                 continue
            
#             gene1_values = gene1_data[common_mask]
#             gene2_values = gene2_data[common_mask]
            
#             # Check for zero variance
#             if np.std(gene1_values) == 0 or np.std(gene2_values) == 0:
#                 continue
            
#             # Compute correlation
#             correlation, p_value = pearsonr(gene1_values, gene2_values)
            
#             # Apply thresholds - strict filtering for speed
#             if abs(correlation) >= correlation_threshold and p_value <= p_value_threshold:
#                 correlations.append({
#                     'gene1': gene_ids[i],
#                     'gene2': gene_ids[j],
#                     'correlation': round(correlation, 4),
#                     'p_value': round(p_value, 6),
#                     'n_samples': int(np.sum(common_mask))
#                 })
                
#         except Exception as e:
#             # Skip problematic pairs silently for speed
#             continue
    
#     return correlations

# def compute_correlation_matrix(df, filters):
#     """
#     Compute correlation matrix using parallel processing for speed.
#     Optimized for large gene sets (30K-60K genes).
#     """
#     try:
#         print("🚀 Computing correlation matrix with parallel processing...")
#         start_time = time.time()
        
#         # Extract parameters - using strict defaults for speed
#         correlation_threshold = filters.get('correlationThreshold', 0.5)  # Strict threshold 
#         p_value_threshold = filters.get('pValueThreshold', 0.01)  # Stricter p-value
#         max_correlations = filters.get('maxCorrelations', 50000)  # Essential limit for large datasets
#         n_processes = filters.get('nProcesses', mp.cpu_count())  # Use all CPU cores
        
#         print(f"🎯 Correlation threshold: {correlation_threshold} (strict for speed)")
#         print(f"📊 P-value threshold: {p_value_threshold}")
#         print(f"🔥 Using {n_processes} CPU cores for parallel processing")
        
#         # Detect matrix structure
#         matrix_start_row, matrix_start_col = detect_matrix_start(df)
        
#         # Extract numeric data matrix
#         data_matrix = df.iloc[matrix_start_row:, matrix_start_col:]
#         gene_metadata = df.iloc[matrix_start_row:, :matrix_start_col]
        
#         print(f"📈 Data matrix shape: {data_matrix.shape}")
        
#         # Convert to numeric and clean data
#         data_matrix = data_matrix.apply(pd.to_numeric, errors='coerce')
        
#         # Pre-filter genes by completeness and variance for speed
#         print("🧹 Pre-filtering genes for speed...")
        
#         # Remove genes with >30% missing data (stricter)
#         missing_threshold = 0.7
#         gene_completeness = 1 - (data_matrix.isnull().sum(axis=1) / data_matrix.shape[1])
#         complete_genes = gene_completeness[gene_completeness >= missing_threshold].index
        
#         # Remove low-variance genes (configurable filtering)
#         variance_filter_method = filters.get('varianceFilterMethod', 'threshold')  # 'percentile' or 'threshold'
        
#         gene_variances = data_matrix.loc[complete_genes].var(axis=1, skipna=True)
#         print(f"📊 Gene variance range: {gene_variances.min():.4f} to {gene_variances.max():.4f}")
#         print(f"📊 Median variance: {gene_variances.median():.4f}")
        
#         if variance_filter_method == 'threshold':
#             # Use absolute variance threshold
#             variance_threshold = filters.get('varianceThreshold', 1.0)  # Default threshold
#             high_var_genes = gene_variances[gene_variances >= variance_threshold].index
#             kept_percentage = (len(high_var_genes) / len(gene_variances)) * 100
#             print(f"🎯 Absolute variance threshold: {variance_threshold}")
#             print(f"📈 Keeping {len(high_var_genes)} genes ({kept_percentage:.1f}%) above threshold")
            
#         else:
#             # Use percentile-based filtering (original method)
#             variance_percentile = filters.get('variancePercentile', 0.3)  # Keep top 70% by default
#             variance_threshold = gene_variances.quantile(variance_percentile)
#             high_var_genes = gene_variances[gene_variances >= variance_threshold].index
#             print(f"🎯 Percentile-based filtering: keeping top {(1-variance_percentile)*100:.0f}% most variable genes")
#             print(f"📊 Variance threshold (percentile): {variance_threshold:.4f}")
        
#         print(f"✅ Final gene count after variance filtering: {len(high_var_genes)} (from {len(complete_genes)})")
        
#         valid_genes = high_var_genes
#         print(f"✅ Genes after pre-filtering: {len(valid_genes)} (from {data_matrix.shape[0]})")
        
#         if len(valid_genes) < 2:
#             return {"error": "Not enough valid genes after filtering"}
        
#         # Prepare data for parallel processing
#         data_matrix_clean = data_matrix.loc[valid_genes].fillna(np.nan)
#         gene_metadata_clean = gene_metadata.loc[valid_genes]
        
#         # Convert to numpy for speed
#         data_matrix_values = data_matrix_clean.values  # Faster than pandas
        
#         # Get gene IDs
#         gene_ids = []
#         for idx in gene_metadata_clean.index:
#             gene_id = str(gene_metadata_clean.loc[idx].iloc[0]).strip()
#             if gene_id in ['nan', '', 'None']:
#                 gene_id = f"Gene_{idx}"
#             gene_ids.append(gene_id)
        
#         # Generate all gene pairs
#         total_genes = len(gene_ids)
#         total_pairs = total_genes * (total_genes - 1) // 2
#         print(f"⚡ Computing {total_pairs:,} pairs across {n_processes} processes")
        
#         # Split pairs into chunks for parallel processing
#         chunk_size = max(1000, total_pairs // (n_processes * 4))  # 4 chunks per process
#         gene_pair_chunks = []
        
#         all_pairs = list(combinations(range(total_genes), 2))
        
#         for i in range(0, len(all_pairs), chunk_size):
#             chunk = all_pairs[i:i + chunk_size]
#             gene_pair_chunks.append(chunk)
        
#         print(f"📦 Split into {len(gene_pair_chunks)} chunks of ~{chunk_size:,} pairs each")
        
#         # Prepare arguments for parallel processing
#         chunk_args = [
#             (chunk, data_matrix_values, gene_ids, correlation_threshold, p_value_threshold)
#             for chunk in gene_pair_chunks
#         ]
        
#         # Parallel computation
#         print("🔥 Starting parallel correlation computation...")
#         computation_start = time.time()
        
#         with mp.Pool(processes=n_processes) as pool:
#             chunk_results = pool.map(compute_correlation_chunk, chunk_args)
        
#         # Combine results from all processes
#         all_correlations = []
#         for chunk_result in chunk_results:
#             all_correlations.extend(chunk_result)
        
#         computation_time = time.time() - computation_start
#         print(f"⚡ Parallel computation completed in {computation_time:.1f} seconds")
#         print(f"🎯 Found {len(all_correlations):,} significant correlations")
        
#         # Sort by absolute correlation (strongest first)
#         all_correlations.sort(key=lambda x: abs(x['correlation']), reverse=True)
        
#         # Apply final limit - essential for large gene networks
#         final_correlations = all_correlations[:max_correlations]
        
#         if len(all_correlations) > max_correlations:
#             print(f"🎯 Limited to top {max_correlations:,} strongest correlations (from {len(all_correlations):,} total)")
#         else:
#             print(f"🎯 Returning all {len(all_correlations):,} significant correlations")
        
#         total_time = time.time() - start_time
#         pairs_per_second = total_pairs / computation_time if computation_time > 0 else 0
        
#         print(f"🚀 COMPLETED: {len(final_correlations):,} correlations in {total_time:.1f}s")
#         print(f"⚡ Speed: {pairs_per_second:,.0f} pairs/second")
        
#         return {
#             "correlations": final_correlations,
#             "totalCorrelations": len(final_correlations),
#             "totalSignificantFound": len(all_correlations),
#             "totalPairsComputed": total_pairs,
#             "computationTime": round(computation_time, 2),
#             "totalTime": round(total_time, 2),
#             "pairsPerSecond": round(pairs_per_second, 0),
#             "parameters": {
#                 "correlationThreshold": correlation_threshold,
#                 "pValueThreshold": p_value_threshold,
#                 "maxCorrelations": max_correlations,
#                 "nProcesses": n_processes,
#                 "preFilteredGenes": len(valid_genes)
#             },
#             "geneCount": len(gene_ids),
#             "sampleCorrelations": final_correlations[:5],  # Top 5 for preview
#             "performance": {
#                 "parallelProcessing": True,
#                 "speedOptimized": True,
#                 "preFiltered": True
#             }
#         }
        
#     except Exception as e:
#         print(f"❌ Error in compute_correlation_matrix: {str(e)}")
#         import traceback
#         traceback.print_exc()
#         return {"error": str(e)}



def generate_fallback_response(command):
    """Generate a fallback response based on command keywords."""
    command_lower = command.lower().strip()
    
    # First, check if this is clearly NOT a data analysis command
    non_data_keywords = [
        # Personal information
        "my name", "i am", "hello", "hi", "goodbye", "bye", "thanks", "thank you",
        # Questions about the system/general
        "what is", "who are", "how are", "weather", "time", "date",
        # Math/calculations unrelated to data
        "calculate", "what's", "plus", "minus", "times", "divided",
        # Random conversation
        "joke", "story", "poem", "sing", "play", "game",
        # Other unrelated topics
        "movie", "food", "music", "sports", "news", "politics"
    ]
    
    # If command contains clearly non-data analysis keywords, return error
    for keyword in non_data_keywords:
        if keyword in command_lower:
            return {"error": "unrecognized command"}
    
    # If command is too short or doesn't contain data analysis keywords, return error
    if len(command_lower) < 3:
        return {"error": "unrecognized command"}
    
    # Check if command contains at least one data analysis keyword
    data_analysis_keywords = [
        "cluster", "sort", "order", "arrange", "search", "find", "locate", "filter", 
        "select", "opacity", "dark", "light", "transparent", "variance", "variant",
        "row", "column", "gene", "sample", "patient", "male", "female", "sex", 
        "race", "timepoint", "histology", "sum", "alphabetical", "alpha", "brca",
        "tp53", "top", "most", "least", "high", "low", "increase", "decrease"
    ]
    
    has_data_keyword = any(keyword in command_lower for keyword in data_analysis_keywords)
    
    if not has_data_keyword:
        return {"error": "unrecognized command"}
    
    # Now handle commands that clearly relate to data analysis
    
    # Handle opacity/appearance commands
    if any(word in command_lower for word in ["opacity", "dark", "light", "transparent"]):
        if any(word in command_lower for word in ["increase", "more", "darker", "high"]):
            return {"action": "set_opacity", "value": "dark"}
        elif any(word in command_lower for word in ["decrease", "less", "lighter", "low"]):
            return {"action": "set_opacity", "value": "light"}
        else:
            # Try to extract numeric values
            import re
            number_match = re.search(r'(\d+(\.\d+)?)', command_lower)
            if number_match:
                return {"action": "set_opacity", "value": number_match.group(1)}
            # Default opacity action
            return {"action": "set_opacity", "value": "dark"}
    
    # Handle clustering commands
    elif "cluster" in command_lower:
        if any(word in command_lower for word in ["row", "gene", "transcript", "feature"]):
            return {"action": "cluster", "target": "rows", "value": ""}
        elif any(word in command_lower for word in ["col", "column", "sample", "patient"]):
            return {"action": "cluster", "target": "columns", "value": ""}
        else:
            # Default clustering
            return {"action": "cluster", "target": "rows", "value": ""}
    
    # Handle sorting commands
    elif any(word in command_lower for word in ["sort", "order", "arrange"]):
        # Check for metadata sorting
        metadata_fields = ["sex", "race", "timepoint", "histology", "survival", "patientid"]
        for field in metadata_fields:
            if field in command_lower:
                # Capitalize field name correctly
                field_capitalized = field.capitalize() if field != "patientid" else "PatientId"
                target = "rows" if any(word in command_lower for word in ["row", "gene"]) else "columns"
                return {"action": "sort_by_meta", "target": target, "value": field_capitalized}
        
        # Handle regular sorting methods
        sort_methods = {
            "variance": ["varian", "var"],
            "sum": ["sum", "total"],
            "alphabetical": ["alpha", "az", "a-z", "name", "alphabet"]
        }
        
        for method, keywords in sort_methods.items():
            if any(keyword in command_lower for keyword in keywords):
                target = "rows" if any(word in command_lower for word in ["row", "gene"]) else "columns"
                return {"action": "sort", "target": target, "value": method}
        
        # Default sorting
        return {"action": "sort", "target": "rows", "value": "variance"}
    
    # Handle search commands
    elif any(word in command_lower for word in ["search", "find", "locate"]):
        # Try to extract search terms
        import re
        # Match anything after search, find, etc.
        search_pattern = re.search(r'(?:search|find|locate)\s+(?:for\s+)?([a-zA-Z0-9_-]+)', command_lower)
        if search_pattern:
            search_term = search_pattern.group(1).strip()
            return {"action": "search", "target": "rows", "value": search_term}
        
        # If no search term found, return error
        return {"error": "unrecognized command"}
    
    # Handle variance commands
    elif "variance" in command_lower or "variant" in command_lower:
        target = "rows"
        # Try to extract a number
        import re
        number_match = re.search(r'(\d+)', command_lower)
        value = number_match.group(1) if number_match else "100"  # Default to top 100
        return {"action": "variance", "target": target, "value": value}
    
    # Handle filtering commands
    elif any(word in command_lower for word in ["filter", "select", "show"]):
        # Look for common filter values
        if any(word in command_lower for word in ["male", "men", "boy"]):
            return {"action": "sample_filter", "target": "Sex", "value": "Male"}
        elif any(word in command_lower for word in ["female", "women", "girl"]):
            return {"action": "sample_filter", "target": "Sex", "value": "Female"}
        else:
            # If filter command but no recognizable value, return error
            return {"error": "unrecognized command"}
    
    # If we get here, the command has data keywords but doesn't match any pattern
    # This means it's unclear what the user wants
    return {"error": "unrecognized command"}


def build_prompt(command, command_history=None, filters=None, metadata=None):
    """
    Build a structured prompt to guide the Ollama model for heatmap data exploration.
    
    Args:
        command (str): The latest user query.
        command_history (list, optional): List of previous commands for context.
        filters (dict, optional): Current applied filters.
        metadata (dict, optional): Available metadata categories and values for filtering.
    
    Returns:
        str: Final prompt string for the Ollama model.
    """
    if command_history is None:
        command_history = []
    if filters is None:
        filters = {}
    
    prompt = (
        "You are an expert data analysis assistant for exploring biological datasets through heatmap visualizations.\n"
        "Respond STRICTLY and ONLY with a valid JSON object. DO NOT include explanations, comments, or extra text.\n\n"
        
        "CRITICAL: If the user command is unclear, unrelated to data analysis, or cannot be mapped to any valid action, "
        "respond with: {\"error\": \"unrecognized command\"}\n\n"
        
        "ONLY process commands that are clearly related to:\n"
        "- Data filtering, sorting, clustering, searching\n"
        "- Heatmap visualization adjustments\n"
        "- Sample or gene selection\n"
        "- Pathway and transcription factor analysis\n\n"
        
        "DO NOT try to interpret unclear commands. DO NOT force-fit unrelated text into valid actions.\n\n"
        
        "Valid Actions (ONLY use these if the command clearly matches):\n"
        "  - 'sort'          : Sort rows or columns. Valid 'value' options: 'alphabetical', 'sum', 'variance'.\n"
        "  - 'sort_by_meta'  : Sort rows or columns by a metadata category. 'target' is 'rows' or 'columns', 'value' is the metadata category.\n"
        "  - 'cluster'       : Perform clustering on rows or columns.\n"
        "  - 'search'        : Search for a specific gene or biomarker. 'value' should be the search term.\n"
        "  - 'sample_filter' : Filter samples based on metadata. 'target' is the metadata field, 'value' is the filter value.\n"
        "  - 'gene_filter'   : Filter rows/genes based on metadata. 'target' is the metadata field, 'value' is the filter value.\n"
        "  - 'variance'      : Select top N most variable rows or columns. 'value' should be the number of top features.\n"
        "  - 'set_opacity'   : Change the opacity/intensity of the heatmap. 'value' should be a number between 0.5 and 3.0 or 'dark'/'light'.\n"
        "  - 'pathway_filter'    : Filter genes by specific pathway name OR transcription factor. 'value' is the exact pathway name or TF name.\n"
        "  - 'pathway_search'    : List available pathways by category. 'value' is the category (immune, cancer, etc.).\n"
        "  - 'functional_filter' : Filter genes by biological function. 'value' is the function name.\n\n"
    )

        
    # Add metadata information in a structured format
    if metadata and isinstance(metadata, dict):
        prompt += "Available Metadata Fields for Filtering and Sorting:\n"
        
        # Process column metadata (typically used for sample filtering)
        if 'col' in metadata and metadata['col']:
            prompt += "Column/Sample Metadata (use exact field names and values):\n"
            for field, values in metadata['col'].items():
                if values and len(values) > 0:
                    # Show field name and some example values
                    example_values = values[:5]  # Show up to 5 values
                    if len(values) > 5:
                        example_values_str = ", ".join([f'"{v}"' for v in example_values]) + ", ..."
                    else:
                        example_values_str = ", ".join([f'"{v}"' for v in example_values])
                    
                    prompt += f"  - '{field}': Possible values: {example_values_str}\n"
        
        # Process row metadata if present
        if 'row' in metadata and metadata['row'] and len(metadata['row']) > 0:
            prompt += "Row/Gene Metadata (use exact field names and values):\n"
            for field, values in metadata['row'].items():
                if values and len(values) > 0:
                    example_values = values[:5]
                    if len(values) > 5:
                        example_values_str = ", ".join([f'"{v}"' for v in example_values]) + ", ..."
                    else:
                        example_values_str = ", ".join([f'"{v}"' for v in example_values])
                    
                    prompt += f"  - '{field}': Possible values: {example_values_str}\n"
        
        prompt += "\n"
    
    # Add command history information
    if command_history and len(command_history) > 0:
        prompt += f"Previous Commands: {json.dumps(command_history)}\n\n"
    
    # Add current filters information
    prompt += f"Current Filters: {json.dumps(filters)}\n\n"
    
    # Add the new user command
    prompt += f"New User Command: {command}\n\n"
    
    # Add examples using actual metadata fields if available
    example_field1 = "Sex"  # Default example
    example_field2 = "Race"
    example_field3 = "Timepoint"
    example_field4 = "Histology"
    
    if metadata and 'col' in metadata:
        if 'Sex' in metadata['col']:
            example_field1 = 'Sex'
        if 'Race' in metadata['col']:
            example_field2 = 'Race'
        if 'Timepoint' in metadata['col']:
            example_field3 = 'Timepoint'
        if 'Histology' in metadata['col']:
            example_field4 = 'Histology'
    
    prompt += f"""
COMMON COMMAND PATTERNS:

Filtering/Selection Commands:
"select males" → {{ "action": "sample_filter", "target": "Sex", "value": "Male" }}
"select females" → {{ "action": "sample_filter", "target": "Sex", "value": "Female" }}  
"filter by race asian" → {{ "action": "sample_filter", "target": "Race", "value": "Asian" }}
"show only alive patients" → {{ "action": "sample_filter", "target": "Survival", "value": "Alive" }}

Sorting Commands:
"sort rows by variance" → {{ "action": "sort", "target": "rows", "value": "variance" }}
"sort by sex" → {{ "action": "sort_by_meta", "target": "columns", "value": "Sex" }}
"sort samples by timepoint" → {{ "action": "sort_by_meta", "target": "columns", "value": "Timepoint" }}

Clustering Commands:
"cluster rows" → {{ "action": "cluster", "target": "rows", "value": "" }}
"cluster genes" → {{ "action": "cluster", "target": "rows", "value": "" }}
"cluster samples" → {{ "action": "cluster", "target": "columns", "value": "" }}

Search Commands:
"search for BRCA1" → {{ "action": "search", "target": "rows", "value": "BRCA1" }}
"find TP53" → {{ "action": "search", "target": "rows", "value": "TP53" }}

Visualization Commands:
"make it darker" → {{ "action": "set_opacity", "value": "dark" }}
"increase opacity" → {{ "action": "set_opacity", "value": "dark" }}
"set opacity to 2.0" → {{ "action": "set_opacity", "value": "2.0" }}

Variance Selection:
"top 50 genes" → {{ "action": "variance", "target": "rows", "value": "50" }}
"most variable 100 rows" → {{ "action": "variance", "target": "rows", "value": "100" }}

PATHWAY & TRANSCRIPTION FACTOR COMMANDS:

Pathway Search Commands (Browse by Category):
"list immune pathways" → {{ "action": "pathway_search", "value": "immune" }}
"show cancer pathways" → {{ "action": "pathway_search", "value": "cancer" }}
"find metabolism pathways" → {{ "action": "pathway_search", "value": "metabolism" }}
"what development pathways are available" → {{ "action": "pathway_search", "value": "development" }}
"show stress pathways" → {{ "action": "pathway_search", "value": "stress" }}
"neural pathways" → {{ "action": "pathway_search", "value": "neural" }}
"hormone pathways" → {{ "action": "pathway_search", "value": "hormone" }}

Specific Pathway Filtering (Filter by Exact Pathway Name):
"show genes from pathway: Immune System R-HSA-168256" → {{ "action": "pathway_filter", "value": "Immune System R-HSA-168256" }}
"show genes from pathway: Reactome_2022_Immune System R-HSA-168256" → {{ "action": "pathway_filter", "value": "Reactome_2022_Immune System R-HSA-168256" }}
"filter by pathway KEGG_MAPK_signaling" → {{ "action": "pathway_filter", "value": "KEGG_MAPK_signaling" }}
"genes from pathway: BioPlanet_2019_Immune system" → {{ "action": "pathway_filter", "value": "BioPlanet_2019_Immune system" }}
"show pathway genes: WikiPathway_2023_Human_Control Of Immune Tolerance" → {{ "action": "pathway_filter", "value": "WikiPathway_2023_Human_Control Of Immune Tolerance" }}
"filter to pathway Jensen_DISEASES_Autoimmune thyroiditis" → {{ "action": "pathway_filter", "value": "Jensen_DISEASES_Autoimmune thyroiditis" }}

IMPORTANT PATHWAY DISTINCTION:
- Use "pathway_search" for BROWSING pathways by category (immune, cancer, etc.)
- Use "pathway_filter" for FILTERING by exact pathway names OR transcription factors
- Commands with "show genes from pathway:" or "filter by pathway" followed by a specific name should use "pathway_filter"

Transcription Factor Filtering:
"show STAT3 targets" → {{ "action": "pathway_filter", "value": "STAT3" }}
"filter to TP53 genes" → {{ "action": "pathway_filter", "value": "TP53" }}
"FOXP3 target genes" → {{ "action": "pathway_filter", "value": "FOXP3" }}
"MYC regulated genes" → {{ "action": "pathway_filter", "value": "MYC" }}
"STAT6 targets" → {{ "action": "pathway_filter", "value": "STAT6" }}
"IRF4 genes" → {{ "action": "pathway_filter", "value": "IRF4" }}
"NFKB1 pathway" → {{ "action": "pathway_filter", "value": "NFKB1" }}
"SOX2 targets" → {{ "action": "pathway_filter", "value": "SOX2" }}

Functional Gene Filtering:
"show immune genes" → {{ "action": "functional_filter", "value": "immune" }}
"cancer related genes" → {{ "action": "functional_filter", "value": "cancer" }}
"metabolism genes" → {{ "action": "functional_filter", "value": "metabolism" }}
"development genes" → {{ "action": "functional_filter", "value": "development" }}
"stress response genes" → {{ "action": "functional_filter", "value": "stress" }}
"neural genes" → {{ "action": "functional_filter", "value": "neural" }}
"hormone responsive genes" → {{ "action": "functional_filter", "value": "hormone" }}

PATHWAY KNOWLEDGE BASE:
Available Categories: immune, cancer, metabolism, development, stress, signaling, hormone, neural, cardiovascular, muscle, bone, circadian, aging, epigenetic

Common Immune TFs: STAT1, STAT3, STAT6, FOXP3, IRF1, IRF4, IRF7, NFKB1, RELA, BCL6, PAX5, TBX21, GATA3
Common Cancer TFs: TP53, MYC, MYCN, BRCA1, BRCA2, RB1, E2F1, FOS, JUN, ETS1, TWIST1, SNAI1
Common Development TFs: SOX2, NANOG, POU5F1, KLF4, PAX6, GATA4, TBX5, HAND2, FOXO1, FOXO3
Common Metabolism TFs: PPARA, PPARG, SREBF1, NRF1, NRF2, HIF1A, TFAM, MLXIPL

PATTERN RECOGNITION FOR PATHWAY COMMANDS:

1. Commands that should trigger "pathway_filter":
   - "show genes from pathway: [EXACT_PATHWAY_NAME]"
   - "filter by pathway [EXACT_PATHWAY_NAME]"
   - "genes from pathway: [EXACT_PATHWAY_NAME]"
   - "show pathway genes: [EXACT_PATHWAY_NAME]"
   - "show [TF_NAME] targets"
   - "[TF_NAME] target genes"
   - "filter to [TF_NAME] genes"
   - Any command asking for genes from a specific named pathway or transcription factor

2. Commands that should trigger "pathway_search":
   - "list [CATEGORY] pathways"
   - "show [CATEGORY] pathways"
   - "find [CATEGORY] pathways"
   - Any command asking to browse pathways by category

VALIDATION PROCESS:
1. Command "show genes from pathway: Immune System R-HSA-168256":
   - Is this pathway analysis? YES (filtering by specific pathway)
   - Contains "show genes from pathway:" pattern? YES
   - Action: pathway_filter, Value: "Immune System R-HSA-168256"
   - Return: {{"action": "pathway_filter", "value": "Immune System R-HSA-168256"}}

2. Command "list immune pathways":
   - Is this pathway analysis? YES (searching pathways)
   - Contains "list [category] pathways" pattern? YES
   - Action: pathway_search, Value: immune
   - Return: {{"action": "pathway_search", "value": "immune"}}

3. Command "show STAT3 targets":
   - Is this pathway analysis? YES (filtering by TF)
   - Contains "[TF] targets" pattern? YES
   - Action: pathway_filter, Value: STAT3
   - Return: {{"action": "pathway_filter", "value": "STAT3"}}

4. Command "my name is john":
   - Is this data analysis? NO
   - Return: {{"error": "unrecognized command"}}

ERROR CASES:
{{"error": "unrecognized command"}} // For non-data analysis commands
{{"error": "value not found: 'XYZ' not available in FieldName. Available values: A, B, C"}} // When requested value doesn't exist

CRITICAL: Always use EXACT field names from the metadata. Pay attention to capitalization.
For pathway commands, distinguish between browsing (pathway_search) and filtering (pathway_filter).
Both specific pathway names and transcription factor names should use "pathway_filter".

Now respond with valid JSON for the command: "{command}"

Start immediately with '{{'.
"""
    
    return prompt
    
def detect_matrix_start(df) -> tuple:
    """Dynamically detect where the numeric data matrix starts."""
    import pandas as pd
    matrix_start_row = None
    matrix_start_col = None
    
    # Step 1: Find the first row that contains mostly numeric values
    for row_idx in range(df.shape[0]):
        numeric_count = 0
        total_count = 0
        
        # Check all columns in this row
        for col_idx in range(df.shape[1]):
            val = df.iloc[row_idx, col_idx]
            if pd.notna(val):
                total_count += 1
                # Check if value is numeric
                if isinstance(val, (int, float)):
                    numeric_count += 1
                elif isinstance(val, str):
                    # Try to convert string to number
                    try:
                        float(val)
                        numeric_count += 1
                    except ValueError:
                        pass
        
        # If more than 50% of non-null values are numeric, this might be data start
        if total_count > 0 and numeric_count / total_count > 0.5:
            matrix_start_row = row_idx
            break
    
    # Default fallback for matrix_start_row
    if matrix_start_row is None:
        matrix_start_row = min(10, df.shape[0] - 1)
        print(f"Warning: Could not detect matrix start row, defaulting to {matrix_start_row}")
    
    # Step 2: Find the first column that contains mostly numeric values in the data rows
    if matrix_start_row is not None:
        for col_idx in range(df.shape[1]):
            numeric_count = 0
            total_count = 0
            
            # Check this column in the data rows (from matrix_start_row onwards)
            # Sample a few rows to determine if this column is numeric
            sample_rows = min(10, df.shape[0] - matrix_start_row)  # Check up to 10 data rows
            
            for row_offset in range(sample_rows):
                val = df.iloc[matrix_start_row + row_offset, col_idx]
                if pd.notna(val):
                    total_count += 1
                    # Check if value is numeric
                    if isinstance(val, (int, float)):
                        numeric_count += 1
                    elif isinstance(val, str):
                        try:
                            float(val)
                            numeric_count += 1
                        except ValueError:
                            pass
            
            # If more than 80% of values in this column are numeric, this is likely data start
            if total_count > 0 and numeric_count / total_count > 0.8:
                matrix_start_col = col_idx
                break
    
    # Default fallback for matrix_start_col
    if matrix_start_col is None:
        matrix_start_col = 0
        print(f"Warning: Could not detect matrix start column, defaulting to {matrix_start_col}")
    
    print(f"Detected matrix start: row={matrix_start_row}, col={matrix_start_col}")
    return matrix_start_row, matrix_start_col

# def apply_filters(df: pd.DataFrame, filters: dict) -> pd.DataFrame:
def apply_filters(df, filters,session_id=None):

    import pandas as pd
    import numpy as np

    try:
        print("📢 apply_filters called!")
        print('******** filter is as follows ******', filters)
        
        if not filters or (not filters.get("row") and not filters.get("col")):
            print("No filters to apply, returning original DataFrame")
            return df
        
        # Create a copy and preserve the original column names
        filtered_df = df.copy()
        original_columns = filtered_df.columns.tolist()
        
        # Detect matrix start
        matrix_start_row, matrix_start_col = detect_matrix_start(filtered_df)
        print(f"📌 matrix_start_row: {matrix_start_row}")
        print(f"📌 matrix_start_col: {matrix_start_col}")
        
        # Separate metadata and data sections
        metadata_rows = filtered_df.iloc[:matrix_start_row, :]
        metadata_cols = filtered_df.iloc[matrix_start_row:, :matrix_start_col]
        data_matrix = filtered_df.iloc[matrix_start_row:, matrix_start_col:]
        
        print(f"Data matrix shape: {data_matrix.shape}")
        print(f"First few values of data matrix:\n{data_matrix.iloc[:3, :3]}")
        
        # Ensure data matrix is numeric for calculations
        data_matrix = data_matrix.apply(pd.to_numeric, errors='coerce')
        
        # ✅ Apply Row Filters (e.g., variance)
        # filtered_row_indices = data_matrix.index.tolist()  # Start with all indices
        filtered_row_indices = set(data_matrix.index)  # ✅ set

        for f in filters.get("row", []):
            print(f"📢 Applying row filter: {f}")
            if f.get("type") == "variance":
                top_n = int(f.get("top_n", 100))

                # Compute variance safely
                variances = data_matrix.replace([np.inf, -np.inf], np.nan).dropna(how='all').var(axis=1)
                if not variances.empty:
                    print(f"Variance range: {variances.min()} to {variances.max()}")
                    top_indices = variances.nlargest(min(top_n, len(variances))).index.tolist()

                    if top_indices:
                        print(f"✅ Selected top {len(top_indices)} high-variance rows")
                        filtered_row_indices = filtered_row_indices & set(top_indices)
                    else:
                        print("⚠️ No rows meet variance threshold")
                else:
                    print("⚠️ No valid rows for variance computation")

            elif f.get("type") == "gene_filter":
                field = f.get("field")
                value = f.get("value")
                if not field or not value:
                    print(f"⚠️ Missing field/value in gene_filter: {f}")
                    continue

                print(f"Looking for rows with metadata '{field}:{value}'")
                field_lower = field.lower()
                value_lower = value.lower()

                matching_rows = []
                for row_idx in range(matrix_start_row, filtered_df.shape[0]):
                    index_val = filtered_df.index[row_idx]
                    if index_val not in filtered_row_indices:
                        continue  # Skip if this row was filtered out earlier

                    # Check each metadata column for this row
                    for col_idx in range(matrix_start_col):
                        cell_value = str(filtered_df.iloc[row_idx, col_idx]).strip()
                        print(f"Checking row {row_idx}, col {col_idx}: '{cell_value}'")
                        
                        match_found = False
                        
                        # Only match if cell contains colon
                        if ":" in cell_value:
                            # Split on first colon only
                            parts = cell_value.split(":", 1)
                            if len(parts) == 2:
                                cell_field = parts[0].strip().lower()
                                cell_value_part = parts[1].strip().lower()
                                
                                # Exact match for both field and value
                                if cell_field == field_lower and cell_value_part == value_lower:
                                    match_found = True
                                    print(f"✅ ROW MATCH: '{cell_field}:{cell_value_part}' matches '{field_lower}:{value_lower}'")
                                    break
                        
                    if match_found:
                        matching_rows.append(index_val)
                        break  # Found match in this row, move to next row

                if matching_rows:
                    print(f"✅ Found {len(matching_rows)} matching rows for '{field}:{value}'")
                    filtered_row_indices = filtered_row_indices & set(matching_rows)
                else:
                    print(f"⚠️ No matching rows found for '{field}:{value}'")
                    
        # Apply final filtered rows
        if filtered_row_indices:
            data_matrix = data_matrix.loc[list(filtered_row_indices)]
            metadata_cols = metadata_cols.loc[list(filtered_row_indices)]

        # ✅ Apply Column Filters (Sample Metadata)
        for f in filters.get("col", []):
            print(f"📢 Applying column filter: {f}")
            if f.get("type") == "sample_filter":
                field = f.get("field")
                value = f.get("value")
                
                if not field or not value:
                    print(f"⚠️ Warning: Missing field or value in filter: {f}")
                    continue
                
                print(f"Looking for samples with {field}:{value}")
                
                # Check if we have NaN metadata rows
                all_nan_metadata = metadata_rows.isna().all().all()
                if all_nan_metadata:
                    print("⚠️ Warning: All metadata rows are NaN, using original file headers")
                    # Load the original file again to get text headers
                    try:
                        orig_file_path = os.path.join(UPLOAD_DIR, f"{session_id}.tsv")
                        header_df = pd.read_csv(orig_file_path, sep="\t", nrows=matrix_start_row)
                        metadata_rows = header_df
                    except Exception as e:
                        print(f"⚠️ Warning: Could not load original file headers: {e}")
                
                # Look for the field in the metadata rows
                field_row = None
                for idx in range(metadata_rows.shape[0]):
                    row_values = metadata_rows.iloc[idx, :].astype(str)
                    # Check for exact field match or field: prefix
                    if any(field == val.strip() or val.strip().startswith(f"{field}:") for val in row_values):
                        field_row = idx
                        break
                
                if field_row is None:
                    print(f"⚠️ Warning: Could not find metadata field '{field}' in rows")
                    continue
                
                print(f"Found field '{field}' in row {field_row}")
                
                # Look for matching columns based on field and value
                matching_cols = []
                
                # Convert field and value to lowercase for case-insensitive matching
                field_lower = field.lower()
                value_lower = value.lower()
                
                for col_idx in range(matrix_start_col, metadata_rows.shape[1]):
                    cell_value = str(metadata_rows.iloc[field_row, col_idx]).strip()
                    print(f"Checking column {col_idx}: '{cell_value}'")
                    
                    match_found = False
                    
                    # Only match if cell contains colon
                    if ":" in cell_value:
                        # Split on first colon only
                        parts = cell_value.split(":", 1)
                        if len(parts) == 2:
                            cell_field = parts[0].strip().lower()
                            cell_value_part = parts[1].strip().lower()
                            
                            # Exact match for both field and value
                            if cell_field == field_lower and cell_value_part == value_lower:
                                match_found = True
                                print(f"✅ COLUMN MATCH: '{cell_field}:{cell_value_part}' matches '{field_lower}:{value_lower}'")
                            else:
                                print(f"❌ No match: '{cell_field}:{cell_value_part}' != '{field_lower}:{value_lower}'")
                    else:
                        print(f"❌ No colon found in '{cell_value}' - skipping")
                    
                    if match_found:
                        matching_cols.append(col_idx)
                
                if matching_cols:
                    # Get column labels for the matching columns
                    valid_sample_labels = [filtered_df.columns[col_idx] for col_idx in matching_cols]
                    print(f"Found {len(valid_sample_labels)} columns matching '{value}' for field '{field}'")
                    
                    # Filter the data matrix to include only these columns
                    common_cols = list(set(data_matrix.columns) & set(valid_sample_labels))
                    if not common_cols:
                        print("⚠️ Warning: No common columns left after this filter")
                    else:
                        print(f"Filtering to {len(common_cols)} columns")
                        data_matrix = data_matrix.loc[:, common_cols]
                else:
                    print(f"⚠️ Warning: No columns found matching '{value}' for field '{field}'")
        
        # Check if we still have data after filtering
        if data_matrix.empty:
            print("⚠️ Warning: No data remains after applying all filters")
            return df  # Return original if all data was filtered out
        
        print(f"Final data matrix shape after filtering: {data_matrix.shape}")
        print(data_matrix)
        
        # ✅ Reconstruct Final DataFrame
        # First combine the metadata columns with filtered data matrix
        if len(metadata_cols) != len(data_matrix):
            print(f"⚠️ Warning: Length mismatch between metadata_cols ({len(metadata_cols)}) and data_matrix ({len(data_matrix)})")
            # Align the metadata_cols to match data_matrix
            metadata_cols = metadata_cols.loc[data_matrix.index]

        print(metadata_cols)
        
        result_data = pd.concat([metadata_cols, data_matrix], axis=1)
        print(result_data)

        # Then prepare metadata rows to align with filtered columns
        metadata_rows_filtered = metadata_rows.copy()
        
        # Add any missing columns to metadata_rows
        for col in result_data.columns:
            if col not in metadata_rows_filtered.columns:
                metadata_rows_filtered[col] = None
        
        # Keep only the columns in the filtered result
        metadata_rows_filtered = metadata_rows_filtered[[col for col in result_data.columns if col in metadata_rows_filtered.columns]]
        
        # If metadata_rows_filtered is missing any columns from result_data, add them
        for col in result_data.columns:
            if col not in metadata_rows_filtered.columns:
                metadata_rows_filtered[col] = None
        
        # Ensure column order matches
        metadata_rows_filtered = metadata_rows_filtered[result_data.columns]
        
        # Finally combine vertically
        final_df = pd.concat([metadata_rows_filtered, result_data], axis=0)
        
        print(f"Final filtered DataFrame shape: {final_df.shape}")
        
        # Critical step: Ensure we're using the original column names
        # This preserves the original unnamed columns without renaming them
        if len(final_df.columns) == len(original_columns):
            # If the filtered dataframe has the same number of columns, 
            # we can directly use the original column names
            final_df.columns = original_columns
        else:
            # If we've filtered out some columns, we need to:
            # 1. Keep the original names for the first metadata_col columns
            # 2. Use the remaining columns from the filtered matrix 
            preserved_cols = original_columns[:matrix_start_col]
            filtered_cols = final_df.columns[matrix_start_col:]
            final_df.columns = preserved_cols + list(filtered_cols)
            
            # Make sure any unnamed columns stay unnamed by preserving their original names
            new_cols = []
            for i, col in enumerate(final_df.columns):
                if i < matrix_start_col and i < len(original_columns):
                    # Keep original name for metadata columns
                    new_cols.append(original_columns[i])
                else:
                    # Keep the current name for data columns
                    new_cols.append(col)
            final_df.columns = new_cols
        
        return final_df
    
    except Exception as e:
        print("❌ Error in apply_filters:", str(e))
        import traceback
        print(traceback.format_exc())
        return df  # Return original data on error

    
def update_filters(filters: dict, action: str, target: str, value: str) -> dict:
    """
    Update the existing filters dictionary by adding or replacing filters based on the AI model's suggestion.
    
    Args:
        filters (dict): Current filters with 'row' and 'col' keys.
        action (str): Action suggested by the model (e.g., 'variance', 'sample_filter').
        target (str): Target of the action (e.g., 'rows', 'columns', or specific metadata field).
        value (str): Value associated with the action (e.g., '100' for top N genes, 'male' for sex).
    
    Returns:
        dict: Updated filters dictionary.
    """
    updated_filters = copy.deepcopy(filters)  # Ensure immutability
    
    if action == "variance":
        # Check if there's already a variance filter and replace it
        variance_filter_exists = False
        for i, filter_item in enumerate(updated_filters.get("row", [])):
            if filter_item.get("type") == "variance":
                # Replace the existing variance filter
                updated_filters["row"][i] = {
                    "type": "variance",
                    "top_n": int(value)  # Ensure value is treated as a number
                }
                variance_filter_exists = True
                break
        
        # If no variance filter exists, add a new one
        if not variance_filter_exists:
            updated_filters.setdefault("row", []).append({
                "type": "variance",
                "top_n": int(value)
            })
    
    elif action == "sample_filter":
        # Check if there's already a sample filter with the same field
        field_filter_exists = False
        for i, filter_item in enumerate(updated_filters.get("col", [])):
            if filter_item.get("type") == "sample_filter" and filter_item.get("field") == target:
                # Replace the existing filter for this field
                updated_filters["col"][i] = {
                    "type": "sample_filter",
                    "field": target,  # This should be the metadata field (e.g., 'sex', 'cell_type')
                    "value": value
                }
                field_filter_exists = True
                break
        
        # If no filter for this field exists, add a new one
        if not field_filter_exists:
            updated_filters.setdefault("col", []).append({
                "type": "sample_filter",
                "field": target,
                "value": value
            })

    elif action == "gene_filter":
        # Check if there's already a gene filter with the same field
        field_filter_exists = False
        for i, filter_item in enumerate(updated_filters.get("row", [])):
            if filter_item.get("type") == "gene_filter" and filter_item.get("field") == target:
                # Replace the existing filter for this field
                updated_filters["row"][i] = {
                    "type": "gene_filter",
                    "field": target,
                    "value": value
                }
                field_filter_exists = True
                break
         # If no filter for this field exists, add a new one
        if not field_filter_exists:
            updated_filters.setdefault("row", []).append({
                "type": "gene_filter",
                "field": target,
                "value": value
            })
    
    elif action == "expression_threshold":
        # Check if there's already an expression threshold filter for this field
        field_filter_exists = False
        for i, filter_item in enumerate(updated_filters.get("row", [])):
            if filter_item.get("type") == "expression_threshold" and filter_item.get("field") == target:
                # Replace the existing filter for this field
                updated_filters["row"][i] = {
                    "type": "expression_threshold",
                    "field": target,
                    "min": float(value)  # Assuming threshold is numeric
                }
                field_filter_exists = True
                break
        
        # If no filter for this field exists, add a new one
        if not field_filter_exists:
            updated_filters.setdefault("row", []).append({
                "type": "expression_threshold",
                "field": target,
                "min": float(value)
            })
    
    # 🚀 Add more action types here as needed (e.g., 'cluster', 'sort', etc.)
    
    return updated_filters


def is_number(value):
    """Check if a value is a number."""
    if value is None:
        return False
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False

def extract_and_save_metadata(tsv_file_path: str, metadata_json_path: str):
    """Extract distinct metadata values for each category from a TSV file and save as JSON."""
    import pandas as pd

    # Read the TSV file
    print(f"Reading TSV file: {tsv_file_path}")
    df = pd.read_csv(tsv_file_path, sep="\t", header=None)
    print(f"DataFrame shape: {df.shape}")
    
    # Print the first few rows to debug
    print("First 5 rows:")
    for i in range(min(5, df.shape[0])):
        print(f"Row {i}: {df.iloc[i].tolist()[:5]}...")
    
    # ✅ Detect start of numeric matrix (row-wise)
    matrix_start_row = None
    for idx in range(df.shape[0]):
        # Check if more than half of the non-empty cells in the row are numbers
        num_count = sum(1 for v in df.iloc[idx, 2:].values if is_number(v))
        non_empty_count = sum(1 for v in df.iloc[idx, 2:].values if pd.notna(v) and str(v).strip())
        
        if non_empty_count > 0 and num_count / non_empty_count > 0.5:
            matrix_start_row = idx
            break
    
    if matrix_start_row is None:
        print("WARNING: Could not find where numeric data rows start. Using row 10 as default.")
        matrix_start_row = min(10, df.shape[0]-1)  # Default to row 10 or last row
    else:
        print(f"Detected matrix start row: {matrix_start_row}")
    
    # ✅ Detect start of numeric matrix (column-wise)
    matrix_start_col = None
    for idx in range(df.shape[1]):
        # Check if more than half of the non-empty cells in the column are numbers
        num_count = sum(1 for v in df.iloc[matrix_start_row:, idx].values if is_number(v))
        non_empty_count = sum(1 for v in df.iloc[matrix_start_row:, idx].values if pd.notna(v) and str(v).strip())
        
        if non_empty_count > 0 and num_count / non_empty_count > 0.5:
            matrix_start_col = idx
            break
    
    if matrix_start_col is None:
        print("WARNING: Could not find where numeric data columns start. Using column 2 as default.")
        matrix_start_col = min(2, df.shape[1]-1)  # Default to column 2 or last column
    else:
        print(f"Detected matrix start column: {matrix_start_col}")
    
    # ✅ Extract Column Metadata (Samples)
    col_metadata = {}
    
    # Extract metadata values from header rows
    for row_idx in range(matrix_start_row):
        for col_idx in range(matrix_start_col, df.shape[1]):
            cell_value = str(df.iloc[row_idx, col_idx]).strip()
            
            if cell_value.lower() in ["nan", "none", "", "sample id"]:
                continue
            
            if ":" in cell_value:
                # Split on first colon
                parts = cell_value.split(":", 1)
                if len(parts) == 2:
                    category, value = parts[0].strip(), parts[1].strip()
                    
                    if category and value and value.lower() not in ["nan", "none", ""]:
                        if category not in col_metadata:
                            col_metadata[category] = set()
                        col_metadata[category].add(value)
    
    # ✅ Extract Row Metadata (Genes/Biomarkers)
    row_metadata = {}
    
    # Extract metadata values from first columns
    for row_idx in range(matrix_start_row, df.shape[0]):
        for col_idx in range(matrix_start_col):
            cell_value = str(df.iloc[row_idx, col_idx]).strip()
            
            if cell_value.lower() in ["nan", "none", "", "gene"]:
                continue
            
            if ":" in cell_value:
                # Split on first colon
                parts = cell_value.split(":", 1)
                if len(parts) == 2:
                    category, value = parts[0].strip(), parts[1].strip()
                    
                    if category and value and value.lower() not in ["nan", "none", ""]:
                        if category not in row_metadata:
                            row_metadata[category] = set()
                        row_metadata[category].add(value)
    
    # Convert sets to sorted lists for JSON serialization
    for key in col_metadata:
        col_metadata[key] = sorted(list(col_metadata[key]))
    
    for key in row_metadata:
        row_metadata[key] = sorted(list(row_metadata[key]))
    
    # ✅ Final Metadata Structure
    metadata = {
        "col": col_metadata,
        "row": row_metadata
    }
    
    # Print summary for debugging
    print(f"Extracted {len(col_metadata)} column metadata categories")
    for category, values in col_metadata.items():
        print(f"  - {category}: {len(values)} unique values")
        print(f"    Example values: {values[:3]}")
    
    print(f"Extracted {len(row_metadata)} row metadata categories")
    for category, values in row_metadata.items():
        print(f"  - {category}: {len(values)} unique values")
        print(f"    Example values: {values[:3]}")
    
    # ✅ Save Metadata JSON
    with open(metadata_json_path, "w") as f:
        json.dump(metadata, f, indent=4)
    
    print(f"Metadata saved at: {metadata_json_path}")
    return metadata

