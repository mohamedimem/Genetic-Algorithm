import os
import wget
from goatools.obo_parser import GODag
from uniport_step1 import UniProtGOExtractor

import pandas as pd
import os
from collections import defaultdict

def save_gene_depths_to_csv(gene_name, go_info_results):
    """
    Save gene depth information to separate CSV files for each depth level.
    Only processes cellular_component GO terms.
    
    Args:
        gene_name (str): Name of the gene
        go_info_results (dict): Results from get_go_info function
    """
    if not go_info_results:
        return
    
    # Filter for cellular_component terms only and group by depth
    depth_go_terms = defaultdict(list)
    
    for go_id, info in go_info_results.items():
        if info is None:
            continue
        if info.get('namespace') == 'cellular_component':
            depth = info.get('depth')
            if depth is not None:
                depth_go_terms[depth].append(go_id)
    
    if not depth_go_terms:
        print(f"No cellular_component GO terms found for {gene_name}")
        return
    
    # Process each depth level
    for depth, go_ids in depth_go_terms.items():
        csv_filename = f"depth_{depth}.csv"
        
        # Check if gene already exists in this depth file
        if os.path.exists(csv_filename):
            try:
                existing_df = pd.read_csv(csv_filename, index_col=0)
                if gene_name in existing_df.index:
                    print(f"Gene {gene_name} already exists in {csv_filename}, skipping...")
                    continue
            except Exception as e:
                print(f"Error reading {csv_filename}: {e}")
                # If file is corrupted, we'll recreate it
                existing_df = pd.DataFrame()
        else:
            existing_df = pd.DataFrame()
        
        # Get all unique GO IDs for this depth across all genes
        all_go_ids = set()
        if not existing_df.empty:
            all_go_ids = set(existing_df.columns)
        all_go_ids.update(go_ids)
        all_go_ids = sorted(list(all_go_ids))  # Sort for consistency
        
        # Create binary vector for current gene
        gene_vector = {}
        for go_id in all_go_ids:
            gene_vector[go_id] = 1 if go_id in go_ids else 0
        
        # Add new gene row
        new_row = pd.DataFrame([gene_vector], index=[gene_name])
        
        # Handle padding for existing genes if new GO IDs were added
        if not existing_df.empty:
            # Add missing columns to existing data
            for go_id in all_go_ids:
                if go_id not in existing_df.columns:
                    existing_df[go_id] = 0
            
            # Reorder columns to match new order
            existing_df = existing_df[all_go_ids]
            
            # Combine existing data with new gene
            updated_df = pd.concat([existing_df, new_row])
        else:
            updated_df = new_row
        
        # Save updated DataFrame
        updated_df.to_csv(csv_filename)
        print(f"Updated {csv_filename} with gene {gene_name}")
    
    # Update count_depth.csv
    update_count_depth_csv(gene_name, depth_go_terms)

def update_count_depth_csv(gene_name, depth_go_terms):
    """
    Update the count_depth.csv file with gene depth counts.
    
    Args:
        gene_name (str): Name of the gene
        depth_go_terms (dict): Dictionary with depth as key and list of GO IDs as values
    """
    count_filename = "count_depth.csv"
    
    # Check if gene already exists
    if os.path.exists(count_filename):
        try:
            count_df = pd.read_csv(count_filename, index_col=0)
            if gene_name in count_df.index:
                print(f"Gene {gene_name} already exists in {count_filename}, skipping...")
                return
        except Exception as e:
            print(f"Error reading {count_filename}: {e}")
            count_df = pd.DataFrame()
    else:
        count_df = pd.DataFrame()
    
    # Calculate depth summary for the gene
    depth_summary = {}
    for depth, go_ids in depth_go_terms.items():
        depth_summary[f"depth_{depth}_count"] = len(go_ids)
    
    # Add total count
    depth_summary["total_go_terms"] = sum(len(go_ids) for go_ids in depth_go_terms.values())
    
    # Create new row
    new_row = pd.DataFrame([depth_summary], index=[gene_name])
    
    # Handle padding - add missing columns to existing data
    if not count_df.empty:
        all_columns = set(count_df.columns) | set(depth_summary.keys())
        
        # Add missing columns to existing data
        for col in all_columns:
            if col not in count_df.columns:
                count_df[col] = 0
            if col not in depth_summary:
                new_row[col] = 0
        
        # Reorder columns
        all_columns = sorted(list(all_columns))
        count_df = count_df[all_columns]
        new_row = new_row[all_columns]
        
        # Combine
        updated_count_df = pd.concat([count_df, new_row])
    else:
        updated_count_df = new_row
    
    # Save updated DataFrame
    updated_count_df.to_csv(count_filename)
    print(f"Updated {count_filename} with gene {gene_name}")

def update_gene_go_matrix(gene_name, go_info_results):
    """
    Alternative method name for save_gene_depths_to_csv for compatibility.
    Does the same thing - processes cellular_component GO terms by depth.
    
    Args:
        gene_name (str): Name of the gene
        go_info_results (dict): Results from get_go_info function
    """
    save_gene_depths_to_csv(gene_name, go_info_results)

# Example usage function to demonstrate how to use with your existing code
def process_gene_list(gene_list, extractor, go_dag):
    """
    Process a list of genes and update all CSV files.
    
    Args:
        gene_list (list): List of gene names
        extractor: UniProtGOExtractor instance
        go_dag: GODag instance
    """
    for gene in gene_list:
        print(f"\n--- Processing gene: {gene} ---")
        try:
            # Get GO terms for the gene
            go_cc_terms = extractor.step1_fetch_go_annotations(gene)
            
            if not go_cc_terms:
                print(f"No GO terms found for {gene}")
                continue
            
            # Get GO information
            target_go_ids = list(go_cc_terms)
            go_info_results = get_go_info(target_go_ids, go_dag)
            
            # Update CSV files
            save_gene_depths_to_csv(gene, go_info_results)
            
        except Exception as e:
            print(f"Error processing gene {gene}: {e}")
            continue

def download_go_obo(url="http://current.geneontology.org/ontology/go-basic.obo", filename="go-basic.obo"):
    """
    Downloads the Gene Ontology OBO file if it doesn't already exist.
    """
    if not os.path.exists(filename):
        print(f"Downloading {filename} from {url}...")
        try:
            wget.download(url, filename)
            print(f"\nSuccessfully downloaded {filename}")
        except Exception as e:
            print(f"\nError downloading file: {e}")
            print("Please ensure you have an active internet connection or download 'go-basic.obo' manually.")
            exit()
    else:
        print(f"{filename} already exists. Skipping download.")

def get_go_info(go_ids, go_dag):
    """
    Retrieves the name, namespace (location), and depth (level) for a list of GO IDs.

    Args:
        go_ids (list): A list of Gene Ontology IDs (e.g., ['GO:0005634', 'GO:0006950']).
        go_dag (GODag): The parsed Gene Ontology DAG object from goatools.

    Returns:
        dict: A dictionary where keys are GO IDs and values are dictionaries
              containing 'name', 'namespace', and 'depth'.
              Returns None for GO IDs not found.
    """
    results = {}
    for go_id in go_ids:
        term = go_dag.get(go_id)
        if term:
            # 'namespace' in GO terms indicates the high-level category (BP, MF, CC)
            # 'depth' indicates the level in the hierarchy (distance from root)
            if term.namespace == "cellular_component":
                results[go_id] = {
                    'name': term.name,
                    'namespace': term.namespace,
                    'depth': term.depth
                }
        else:
            results[go_id] = None # Indicate that the GO ID was not found
            print(f"Warning: GO ID '{go_id}' not found in the OBO file.")
    return results

if __name__ == "__main__":
    # Initialize
    extractor = UniProtGOExtractor()

    # Extract GO:CC terms for a gene
    go_cc_terms = extractor.step1_fetch_go_annotations("POLE3")

    # Results: {'GO:0005634', 'GO:0005737', 'GO:0005739'}
    
    # --- Step 1: Download the GO OBO file ---
    download_go_obo(filename="go-basic.obo")
    go_dag = GODag("go-basic.obo")
    # --- Step 2: Load the OBO file into a G

    # --- Step 3: Process POLE3 example ---
    target_go_ids = list(go_cc_terms)
    go_info_results = get_go_info(target_go_ids, go_dag)
    save_gene_depths_to_csv("gene", go_info_results)
            
