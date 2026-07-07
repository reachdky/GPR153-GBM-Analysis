#!/usr/bin/env python3
"""
Pathway Analysis - GPR153 Correlated Genes
Performs pathway enrichment analysis on genes correlated with GPR153
"""

import os
import sys
import pandas as pd
import requests
import time
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def load_correlated_genes():
    """Load GPR153 correlated genes"""
    corr_path = 'data/results/correlation/gpr153_correlated_genes.csv'
    
    if not os.path.exists(corr_path):
        log_message("ERROR: Correlation results not found. Run script 09 first.")
        return None
    
    log_message("Loading GPR153 correlated genes...")
    corr_data = pd.read_csv(corr_path)
    
    log_message(f"Loaded {len(corr_data)} correlated genes")
    return corr_data

def enrichr_analysis(gene_list, gene_set_library='KEGG_2021_Human'):
    """Perform enrichment analysis using Enrichr API"""
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'
    
    # Submit gene list
    genes_str = '\n'.join(gene_list)
    payload = {
        'list': (None, genes_str),
        'description': (None, 'GPR153 correlated genes')
    }
    
    try:
        response = requests.post(f'{ENRICHR_URL}/addList', files=payload)
        if response.status_code != 200:
            return None
        
        user_list_id = response.json()['userListId']
        
        # Get enrichment results
        query_string = f'?userListId={user_list_id}&backgroundType={gene_set_library}'
        response = requests.get(f'{ENRICHR_URL}/enrich{query_string}')
        
        if response.status_code != 200:
            return None
        
        data = response.json()
        
        if gene_set_library not in data:
            return None
        
        return data[gene_set_library]
    
    except Exception as e:
        log_message(f"Enrichr API error: {e}")
        return None

def parse_enrichr_results(results, top_n=20):
    """Parse Enrichr results into DataFrame"""
    if not results:
        return pd.DataFrame()
    
    parsed = []
    for result in results[:top_n]:
        parsed.append({
            'term': result[1],
            'p_value': result[2],
            'adjusted_p_value': result[6],
            'z_score': result[3],
            'combined_score': result[4],
            'genes': ';'.join(result[5]),
            'n_genes': len(result[5])
        })
    
    return pd.DataFrame(parsed)

def main():
    """Main pathway analysis pipeline"""
    log_message("Starting pathway analysis")
    print("=" * 70)
    print("PATHWAY ENRICHMENT ANALYSIS")
    print("=" * 70)
    
    # Load correlated genes
    corr_data = load_correlated_genes()
    if corr_data is None:
        return
    
    # Create output directory
    output_dir = 'data/results/pathways'
    os.makedirs(output_dir, exist_ok=True)
    
    # Get gene lists
    all_genes = corr_data['gene'].tolist()
    positive_corr = corr_data[corr_data['correlation'] > 0]['gene'].tolist()
    negative_corr = corr_data[corr_data['correlation'] < 0]['gene'].tolist()
    
    print(f"\nGene sets:")
    print(f"  All correlated genes: {len(all_genes)}")
    print(f"  Positively correlated: {len(positive_corr)}")
    print(f"  Negatively correlated: {len(negative_corr)}")
    
    # KEGG pathway enrichment
    print("\n" + "=" * 70)
    print("KEGG PATHWAY ENRICHMENT (All genes)")
    print("=" * 70)
    
    log_message("Querying Enrichr API for KEGG pathways...")
    kegg_results = enrichr_analysis(all_genes[:500], 'KEGG_2021_Human')
    time.sleep(1)  # Be nice to API
    
    if kegg_results:
        kegg_df = parse_enrichr_results(kegg_results, top_n=20)
        
        if len(kegg_df) > 0:
            print(f"\nTop 10 KEGG pathways:")
            for idx, row in kegg_df.head(10).iterrows():
                print(f"  {row['term']}")
                print(f"    p-value: {row['p_value']:.2e}, genes: {row['n_genes']}")
            
            kegg_df.to_csv(f'{output_dir}/kegg_pathways.csv', index=False)
            print(f"\n  ✓ Saved: kegg_pathways.csv")
        else:
            print("  No significant pathways found")
    else:
        print("  ✗ API request failed")
    
    # GO Biological Process
    print("\n" + "=" * 70)
    print("GO BIOLOGICAL PROCESS ENRICHMENT")
    print("=" * 70)
    
    log_message("Querying Enrichr API for GO Biological Process...")
    go_bp_results = enrichr_analysis(all_genes[:500], 'GO_Biological_Process_2021')
    time.sleep(1)
    
    if go_bp_results:
        go_bp_df = parse_enrichr_results(go_bp_results, top_n=20)
        
        if len(go_bp_df) > 0:
            print(f"\nTop 10 GO Biological Processes:")
            for idx, row in go_bp_df.head(10).iterrows():
                print(f"  {row['term']}")
                print(f"    p-value: {row['p_value']:.2e}, genes: {row['n_genes']}")
            
            go_bp_df.to_csv(f'{output_dir}/go_biological_process.csv', index=False)
            print(f"\n  ✓ Saved: go_biological_process.csv")
        else:
            print("  No significant processes found")
    else:
        print("  ✗ API request failed")
    
    # GO Molecular Function
    print("\n" + "=" * 70)
    print("GO MOLECULAR FUNCTION ENRICHMENT")
    print("=" * 70)
    
    log_message("Querying Enrichr API for GO Molecular Function...")
    go_mf_results = enrichr_analysis(all_genes[:500], 'GO_Molecular_Function_2021')
    time.sleep(1)
    
    if go_mf_results:
        go_mf_df = parse_enrichr_results(go_mf_results, top_n=20)
        
        if len(go_mf_df) > 0:
            print(f"\nTop 10 GO Molecular Functions:")
            for idx, row in go_mf_df.head(10).iterrows():
                print(f"  {row['term']}")
                print(f"    p-value: {row['p_value']:.2e}, genes: {row['n_genes']}")
            
            go_mf_df.to_csv(f'{output_dir}/go_molecular_function.csv', index=False)
            print(f"\n  ✓ Saved: go_molecular_function.csv")
        else:
            print("  No significant functions found")
    else:
        print("  ✗ API request failed")
    
    # Reactome pathways
    print("\n" + "=" * 70)
    print("REACTOME PATHWAY ENRICHMENT")
    print("=" * 70)
    
    log_message("Querying Enrichr API for Reactome pathways...")
    reactome_results = enrichr_analysis(all_genes[:500], 'Reactome_2022')
    time.sleep(1)
    
    if reactome_results:
        reactome_df = parse_enrichr_results(reactome_results, top_n=20)
        
        if len(reactome_df) > 0:
            print(f"\nTop 10 Reactome pathways:")
            for idx, row in reactome_df.head(10).iterrows():
                print(f"  {row['term']}")
                print(f"    p-value: {row['p_value']:.2e}, genes: {row['n_genes']}")
            
            reactome_df.to_csv(f'{output_dir}/reactome_pathways.csv', index=False)
            print(f"\n  ✓ Saved: reactome_pathways.csv")
        else:
            print("  No significant pathways found")
    else:
        print("  ✗ API request failed")
    
    print("\n" + "=" * 70)
    log_message("Pathway analysis complete")
    print("=" * 70)
    print("\nNote: Pathway enrichment requires internet connection.")
    print("Results depend on Enrichr API availability.")

if __name__ == '__main__':
    main()
