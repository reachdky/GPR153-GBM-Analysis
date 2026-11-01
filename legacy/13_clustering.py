#!/usr/bin/env python3
"""
Clustering Analysis
Performs unsupervised clustering on GBM samples based on gene expression
"""

import os
import sys
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def load_data():
    """Load expression data"""
    data_path = 'data/processed/combined_gbm_normal_expression.csv'
    
    log_message("Loading expression data...")
    data = pd.read_csv(data_path, index_col=0)
    
    # Filter for GBM samples only
    gbm_data = data[data['group'] == 'GBM'].copy()
    gbm_data = gbm_data.drop('group', axis=1)
    
    log_message(f"Loaded GBM data: {gbm_data.shape}")
    return gbm_data

def select_variable_genes(data, n_genes=1000):
    """Select most variable genes"""
    # Calculate variance for each gene
    gene_var = data.var(axis=0)
    
    # Select top variable genes
    top_genes = gene_var.nlargest(n_genes).index
    
    return data[top_genes]

def perform_pca(data, n_components=10):
    """Perform PCA"""
    # Standardize
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    
    # PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data_scaled)
    
    # Create DataFrame
    pca_df = pd.DataFrame(
        pca_result,
        index=data.index,
        columns=[f'PC{i+1}' for i in range(n_components)]
    )
    
    return pca_df, pca.explained_variance_ratio_

def perform_kmeans(data, n_clusters=3):
    """Perform K-means clustering"""
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(data)
    
    return clusters, kmeans

def analyze_cluster_characteristics(original_data, clusters, genes_of_interest):
    """Analyze characteristics of each cluster"""
    results = []
    
    for cluster_id in np.unique(clusters):
        cluster_mask = clusters == cluster_id
        cluster_samples = original_data[cluster_mask]
        
        char = {
            'cluster': cluster_id,
            'n_samples': cluster_mask.sum()
        }
        
        # Mean expression of genes of interest
        for gene in genes_of_interest:
            if gene in original_data.columns:
                char[f'{gene}_mean'] = cluster_samples[gene].mean()
        
        results.append(char)
    
    return pd.DataFrame(results)

def main():
    """Main clustering pipeline"""
    log_message("Starting clustering analysis")
    print("=" * 70)
    print("CLUSTERING ANALYSIS")
    print("=" * 70)
    
    # Load data
    data = load_data()
    
    # Create output directory
    output_dir = 'data/results/ml'
    os.makedirs(output_dir, exist_ok=True)
    
    # Select variable genes
    print("\n" + "=" * 70)
    print("GENE SELECTION")
    print("=" * 70)
    
    var_data = select_variable_genes(data, n_genes=1000)
    print(f"  Selected {var_data.shape[1]} most variable genes")
    
    # PCA
    print("\n" + "=" * 70)
    print("PRINCIPAL COMPONENT ANALYSIS")
    print("=" * 70)
    
    pca_data, variance_ratio = perform_pca(var_data, n_components=10)
    
    print(f"\nExplained variance by component:")
    for i, var in enumerate(variance_ratio, 1):
        print(f"  PC{i}: {var*100:.2f}%")
    
    print(f"\nCumulative variance (PC1-3): {variance_ratio[:3].sum()*100:.2f}%")
    
    # Save PCA results
    pca_data.to_csv(f'{output_dir}/pca_results.csv')
    
    # Save variance
    var_df = pd.DataFrame({
        'component': [f'PC{i+1}' for i in range(len(variance_ratio))],
        'variance_explained': variance_ratio
    })
    var_df.to_csv(f'{output_dir}/pca_variance.csv', index=False)
    
    # K-means clustering
    print("\n" + "=" * 70)
    print("K-MEANS CLUSTERING")
    print("=" * 70)
    
    # Try different numbers of clusters
    for n_clusters in [2, 3, 4]:
        clusters, kmeans = perform_kmeans(pca_data.iloc[:, :3], n_clusters=n_clusters)
        
        print(f"\n{n_clusters} clusters:")
        unique, counts = np.unique(clusters, return_counts=True)
        for cluster_id, count in zip(unique, counts):
            print(f"  Cluster {cluster_id}: {count} samples")
    
    # Use 3 clusters for final analysis
    clusters, kmeans = perform_kmeans(pca_data.iloc[:, :3], n_clusters=3)
    
    # Save cluster assignments
    cluster_df = pd.DataFrame({
        'sample_id': data.index,
        'cluster': clusters
    })
    cluster_df.to_csv(f'{output_dir}/cluster_assignments.csv', index=False)
    
    # Analyze cluster characteristics
    print("\n" + "=" * 70)
    print("CLUSTER CHARACTERISTICS")
    print("=" * 70)
    
    genes_of_interest = ['GPR153', 'MGMT', 'TP53', 'EGFR', 'IDH1', 'PTEN']
    cluster_char = analyze_cluster_characteristics(data, clusters, genes_of_interest)
    
    print("\nMean expression by cluster:")
    print(cluster_char.to_string(index=False))
    
    cluster_char.to_csv(f'{output_dir}/cluster_characteristics.csv', index=False)
    
    print(f"\n✓ Saved: pca_results.csv")
    print(f"✓ Saved: cluster_assignments.csv")
    print(f"✓ Saved: cluster_characteristics.csv")
    
    print("\n" + "=" * 70)
    log_message("Clustering analysis complete")
    print("=" * 70)

if __name__ == '__main__':
    main()
