import pandas as pd
import numpy as np
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
import threading
import os
from typing import Dict, List, Tuple, Optional, Union
import time
import warnings
import json
warnings.filterwarnings('ignore')

class GeneSimilarityAnalyzer:
    def __init__(self, data_dir: str = "./", sigma: float = 1.2, min_frequency: int = 2):
        """
        Initialize the Gene Similarity Analyzer
        
        Args:
            data_dir: Directory containing the CSV files
            sigma: Gaussian kernel parameter
            min_frequency: Minimum frequency for Adamic/Adar to avoid log(1)
        """
        self.data_dir = data_dir
        self.sigma = sigma
        self.min_frequency = min_frequency
        self.depth_data = {}  # Store depth_X.csv data
        self.go_frequencies = {}  # Store GO term frequencies per level
        self.gene_levels = {}  # Track which levels each gene appears in
        self.lock = threading.Lock()
        self.checkpoint_file = "analysis_checkpoint.json"
        self.temp_output_file = "temp_results.csv"
        
    def save_checkpoint(self, completed_chunks: int, total_chunks: int, results_so_far: List):
        """Save progress checkpoint"""
        checkpoint_data = {
            "completed_chunks": completed_chunks,
            "total_chunks": total_chunks,
            "total_completed_pairs": len(results_so_far),
            "timestamp": time.time()
        }
        
        with self.lock:
            # Save checkpoint info
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f)
            
            # Save results so far to temporary file
            if results_so_far:
                temp_df = pd.DataFrame(results_so_far, 
                                     columns=['Gene1', 'Gene2', 'jaccard_max', 
                                             'adamic_max', 'gaussian_max'])
                temp_df.to_csv(self.temp_output_file, sep=';', index=False)
                
            print(f"✓ Checkpoint saved: {completed_chunks}/{total_chunks} chunks complete "
                  f"({len(results_so_far):,} pairs processed)")
    
    def load_checkpoint(self):
        """Load existing checkpoint if available"""
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, 'r') as f:
                    checkpoint_data = json.load(f)
                
                # Load temporary results if they exist
                completed_results = []
                if os.path.exists(self.temp_output_file):
                    temp_df = pd.read_csv(self.temp_output_file, sep=';')
                    completed_results = temp_df.to_dict('records')
                
                return checkpoint_data, completed_results
            except:
                print("Warning: Could not load checkpoint file")
                return None, []
        return None, []
    
    def cleanup_checkpoint(self):
        """Clean up checkpoint files after successful completion"""
        try:
            if os.path.exists(self.checkpoint_file):
                os.remove(self.checkpoint_file)
            if os.path.exists(self.temp_output_file):
                os.remove(self.temp_output_file)
            print("✓ Checkpoint files cleaned up")
        except:
            print("Warning: Could not clean up checkpoint files")
    
    def debug_data_structure(self):
        """Debug function to understand your data structure"""
        print("\n=== DATA STRUCTURE DEBUG ===")
        
        # Check count_depth.csv
        count_file = os.path.join(self.data_dir, "count_depth.csv")
        if os.path.exists(count_file):
            print(f"\n1. ANALYZING: {count_file}")
            count_df = pd.read_csv(count_file)
            print(f"   Shape: {count_df.shape}")
            print(f"   Columns: {count_df.columns.tolist()}")
            print(f"   First 3 rows:")
            print(count_df.head(3))
        
        # Check depth files
        print(f"\n2. ANALYZING DEPTH FILES:")
        for depth in range(1, 11):
            file_path = os.path.join(self.data_dir, f"depth_{depth}.csv")
            if os.path.exists(file_path):
                df = pd.read_csv(file_path)
                print(f"   depth_{depth}.csv: {df.shape} - Columns: {df.columns[:5].tolist()}...")
                break  # Just show one example
        
        print("=== END DEBUG ===\n")
    
    def load_data(self):
        """Load all depth_X.csv files and calculate GO frequencies"""
        print("Loading data files...")
        
        # First, try to load count_depth.csv to get ALL genes (including those with no annotations)
        count_file = os.path.join(self.data_dir, "count_depth.csv")
        if os.path.exists(count_file):
            print(f"Loading {count_file}...")
            count_df = pd.read_csv(count_file)
            
            # Debug: Print column names to understand the structure
            print(f"Columns in count_depth.csv: {count_df.columns.tolist()}")
            
            # Try different possible column names for genes
            gene_column = None
            possible_gene_columns = ['Gene', 'gene', 'GENE', 'gene_name', 'Gene_name', 'symbol']
            
            for col in possible_gene_columns:
                if col in count_df.columns:
                    gene_column = col
                    break
            
            if gene_column:
                all_genes = count_df[gene_column].tolist()
                print(f"Found {len(all_genes)} total genes in count_depth.csv using column '{gene_column}'")
                
                # Initialize all genes with empty level lists
                for gene in all_genes:
                    self.gene_levels[gene] = []
            else:
                # If no recognizable gene column, use the first column or index
                if count_df.index.name or len(count_df.columns) > 0:
                    # Try using index if it looks like gene names
                    if count_df.index.name and count_df.index.name.lower() in ['gene', 'genes']:
                        all_genes = count_df.index.tolist()
                        print(f"Using index as gene names: {len(all_genes)} genes found")
                    else:
                        # Use first column as gene names
                        all_genes = count_df.iloc[:, 0].tolist()
                        print(f"Using first column as gene names: {len(all_genes)} genes found")
                    
                    # Initialize all genes with empty level lists
                    for gene in all_genes:
                        self.gene_levels[gene] = []
                else:
                    print("Warning: Could not identify gene column in count_depth.csv")
                    print("Will only use genes found in depth files.")
        else:
            print("Warning: count_depth.csv not found. Will only use genes from depth files.")
        
        # Load depth files (1-10)
        for depth in range(1, 11):
            file_path = os.path.join(self.data_dir, f"depth_{depth}.csv")
            if os.path.exists(file_path):
                print(f"Loading depth_{depth}.csv...")
                
                try:
                    # Try to read with index_col=0 (gene names as index)
                    df = pd.read_csv(file_path, index_col=0)
                except:
                    # If that fails, read normally and use first column as index
                    df = pd.read_csv(file_path)
                    df = df.set_index(df.columns[0])
                
                self.depth_data[depth] = df
                
                # Calculate GO term frequencies for Adamic/Adar
                self.go_frequencies[depth] = {}
                for go_term in df.columns:
                    freq = df[go_term].sum()
                    self.go_frequencies[depth][go_term] = max(freq, 1)  # Avoid 0
                
                # Track which genes appear at this level
                for gene in df.index:
                    if gene not in self.gene_levels:
                        self.gene_levels[gene] = []
                    self.gene_levels[gene].append(depth)
                    
                print(f"  - Loaded {len(df)} genes, {len(df.columns)} GO terms")
            else:
                print(f"Warning: {file_path} not found")
        
        print(f"Data loading complete. Total unique genes: {len(self.gene_levels)}")
        
        # Report genes with no annotations
        no_annotation_genes = [gene for gene, levels in self.gene_levels.items() if not levels]
        if no_annotation_genes:
            print(f"Warning: {len(no_annotation_genes)} genes have no GO annotations")
            print(f"These will have similarity = None with all other genes")
        
    def jaccard_similarity(self, vec_a: np.array, vec_b: np.array) -> float:
        """Calculate Jaccard similarity between two binary vectors"""
        intersection = np.sum(vec_a & vec_b)
        union = np.sum(vec_a | vec_b)
        return intersection / union if union > 0 else 0.0
    
    def adamic_adar_similarity(self, vec_a: np.array, vec_b: np.array, 
                             go_terms: List[str], level: int) -> float:
        """Calculate Adamic/Adar similarity between two binary vectors"""
        score = 0.0
        frequencies = self.go_frequencies[level]
        
        for i, go_term in enumerate(go_terms):
            if vec_a[i] == 1 and vec_b[i] == 1:  # Both genes have this GO term
                freq = frequencies[go_term]
                if freq >= self.min_frequency:
                    score += 1.0 / math.log(freq)
        
        return score
    
    def gaussian_similarity(self, vec_a: np.array, vec_b: np.array) -> float:
        """Calculate Gaussian kernel similarity between two vectors"""
        euclidean_dist_sq = np.sum((vec_a - vec_b) ** 2)
        return math.exp(-euclidean_dist_sq / (2 * self.sigma ** 2))
    
    def calculate_similarity_at_level(self, gene_a: str, gene_b: str, level: int) -> Tuple[float, float, float]:
        """Calculate all three similarities for a gene pair at a specific level"""
        df = self.depth_data[level]
        
        # Get vectors for both genes
        vec_a = df.loc[gene_a].values.astype(int)
        vec_b = df.loc[gene_b].values.astype(int)
        go_terms = df.columns.tolist()
        
        # Calculate all three similarities
        jaccard = self.jaccard_similarity(vec_a, vec_b)
        adamic = self.adamic_adar_similarity(vec_a, vec_b, go_terms, level)
        gaussian = self.gaussian_similarity(vec_a, vec_b)
        
        return jaccard, adamic, gaussian
    
    def calculate_gene_pair_similarity(self, gene_a: str, gene_b: str) -> Tuple[str, str, Union[float, None], Union[float, None], Union[float, None]]:
        """Calculate maximum similarity across all common levels for a gene pair"""
        # Handle missing genes gracefully - return None instead of 0.0
        if gene_a not in self.gene_levels or gene_b not in self.gene_levels:
            return gene_a, gene_b, None, None, None
        
        # Find common levels where both genes appear
        levels_a = set(self.gene_levels.get(gene_a, []))
        levels_b = set(self.gene_levels.get(gene_b, []))
        common_levels = levels_a & levels_b
        
        # If no common levels, return None instead of 0.0
        if not common_levels:
            return gene_a, gene_b, None, None, None
        
        max_jaccard = 0.0
        max_adamic = 0.0
        max_gaussian = 0.0
        
        # Calculate similarities at each common level and take maximum
        for level in common_levels:
            try:
                jaccard, adamic, gaussian = self.calculate_similarity_at_level(gene_a, gene_b, level)
                max_jaccard = max(max_jaccard, jaccard)
                max_adamic = max(max_adamic, adamic)
                max_gaussian = max(max_gaussian, gaussian)
            except KeyError:
                # Gene not found at this level (shouldn't happen with proper preprocessing)
                continue
        
        return gene_a, gene_b, max_jaccard, max_adamic, max_gaussian
    
    def process_gene_pairs_chunk(self, gene_pairs_chunk: List[Tuple[str, str]], 
                                chunk_id: int) -> List[Tuple[str, str, Union[float, None], Union[float, None], Union[float, None]]]:
        """Process a chunk of gene pairs from the CSV file"""
        results = []
        total_pairs = len(gene_pairs_chunk)
        
        for i, (gene_a, gene_b) in enumerate(gene_pairs_chunk):
            if i % 1000 == 0:  # Progress reporting
                with self.lock:
                    print(f"Chunk {chunk_id}: Processed {i}/{total_pairs} pairs")
            
            result = self.calculate_gene_pair_similarity(gene_a, gene_b)
            results.append(result)
        
        with self.lock:
            print(f"Chunk {chunk_id}: Completed {total_pairs} pairs")
        
        return results
    
    def process_gene_chunk(self, gene_pairs_chunk: List[Tuple[str, str]], 
                          chunk_id: int) -> List[Tuple[str, str, Union[float, None], Union[float, None], Union[float, None]]]:
        """Process a chunk of gene pairs in a separate thread (legacy method for full analysis)"""
        results = []
        total_pairs = len(gene_pairs_chunk)
        
        for i, (gene_a, gene_b) in enumerate(gene_pairs_chunk):
            if i % 1000 == 0:  # Progress reporting
                with self.lock:
                    print(f"Chunk {chunk_id}: Processed {i}/{total_pairs} pairs")
            
            result = self.calculate_gene_pair_similarity(gene_a, gene_b)
            results.append(result)
        
        with self.lock:
            print(f"Chunk {chunk_id}: Completed {total_pairs} pairs")
        
        return results
    
    def load_gene_pairs_from_csv(self, input_file: str) -> pd.DataFrame:
        """Load gene pairs from existing CSV file"""
        print(f"Loading gene pairs from {input_file}...")
        
        # Read the CSV with semicolon separator
        df = pd.read_csv(input_file, sep=';')
        
        print(f"Loaded {len(df)} gene pairs from CSV")
        print(f"Columns found: {df.columns.tolist()}")
        
        # Check if similarity columns already exist
        similarity_columns = ['jaccard_max', 'adamic_max', 'gaussian_max']
        existing_columns = [col for col in similarity_columns if col in df.columns]
        
        if existing_columns:
            print(f"🔍 Found existing similarity columns: {existing_columns}")
            
            # Check how many rows already have values
            has_values_mask = True
            for col in existing_columns:
                # Consider a row "complete" if it has non-null values in all existing similarity columns
                # Note: We now allow 0.0 as a valid value, only skip if completely null
                has_values_mask = has_values_mask & df[col].notna()
            
            completed_rows = has_values_mask.sum()
            total_rows = len(df)
            
            print(f"📊 Analysis status:")
            print(f"   - Total pairs: {total_rows:,}")
            print(f"   - Already calculated: {completed_rows:,}")
            print(f"   - Need calculation: {(total_rows - completed_rows):,}")
            
            if completed_rows > 0:
                skip_choice = input("Skip already calculated pairs? (y/n): ").lower().strip()
                if skip_choice == 'y':
                    print(f"✅ Will skip {completed_rows:,} pairs that already have similarity values")
                    return df, has_values_mask
                else:
                    print("🔄 Will recalculate all pairs")
        
        return df, None
    
    def get_all_gene_pairs(self) -> List[Tuple[str, str]]:
        """Generate all unique gene pairs in consistent alphabetical order"""
        all_genes = sorted(list(self.gene_levels.keys()))  # Sort for consistent ordering
        print(f"Generating gene pairs from {len(all_genes)} genes...")
        
        # Generate ordered pairs: ensure gene_a < gene_b alphabetically
        gene_pairs = [(gene_a, gene_b) for gene_a, gene_b in combinations(all_genes, 2)]
        print(f"Total gene pairs to process: {len(gene_pairs):,}")
        
        return gene_pairs
    
    def run_analysis_from_csv(self, input_csv_file: str, output_file: str = "gene_similarity_results_enhanced.csv", 
                             chunk_size: int = 10000, max_workers: int = None, 
                             save_frequency: int = 5):
        """Run similarity analysis with checkpoint/resume functionality and skip-if-exists"""
        
        # Check for existing checkpoint
        checkpoint_data, completed_results = self.load_checkpoint()
        
        if checkpoint_data:
            print(f"\n🔄 RESUMING FROM CHECKPOINT:")
            print(f"   Previously completed: {checkpoint_data['completed_chunks']}/{checkpoint_data['total_chunks']} chunks")
            print(f"   Previously processed: {checkpoint_data['total_completed_pairs']:,} pairs")
            
            resume_choice = input("Resume from checkpoint? (y/n): ").lower().strip()
            if resume_choice != 'y':
                print("Starting fresh analysis...")
                completed_results = []
                checkpoint_data = None
        
        # Load GO annotation data
        print("Loading GO annotation data...")
        self.load_data()
        
        # Load your gene pairs from CSV and check for existing values
        input_df, skip_mask = self.load_gene_pairs_from_csv(input_csv_file)
        
        if len(input_df) == 0:
            print("No gene pairs found in CSV file!")
            return
        
        # Determine which pairs need processing
        if skip_mask is not None:
            # Only process pairs that don't have values yet
            pairs_to_process = input_df[~skip_mask]
            print(f"🎯 Processing {len(pairs_to_process):,} pairs that need calculation")
        else:
            # Process all pairs
            pairs_to_process = input_df
            print(f"🎯 Processing all {len(pairs_to_process):,} pairs")
        
        if len(pairs_to_process) == 0:
            print("🎉 All pairs already have similarity values! Nothing to calculate.")
            return input_df
        
        # Extract gene pairs that need processing
        gene_pairs_to_calc = [(row['Gene1'], row['Gene2']) for _, row in pairs_to_process.iterrows()]
        
        # Set default number of workers
        if max_workers is None:
            max_workers = min(os.cpu_count(), 8)
        
        # Split gene pairs into chunks
        chunks = [gene_pairs_to_calc[i:i + chunk_size] 
                 for i in range(0, len(gene_pairs_to_calc), chunk_size)]
        
        total_chunks = len(chunks)
        print(f"Created {total_chunks} chunks (chunk size: {chunk_size:,})")
        
        # Determine starting point based on checkpoint
        start_chunk = 0
        all_results = completed_results.copy() if completed_results else []
        
        if checkpoint_data and len(gene_pairs_to_calc) > 0:
            start_chunk = checkpoint_data['completed_chunks']
            chunks = chunks[start_chunk:]  # Skip already completed chunks
            print(f"Resuming from chunk {start_chunk + 1}/{total_chunks}")
        
        print(f"Processing with {max_workers} threads...")
        print(f"Auto-save frequency: Every {save_frequency} chunks")
        
        # Process chunks in parallel with checkpointing
        start_time = time.time()
        completed_chunks = start_chunk
        
        try:
            if len(chunks) > 0:  # Only process if there are chunks to process
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    # Process chunks in batches for checkpointing
                    chunk_batches = [chunks[i:i + save_frequency] 
                                   for i in range(0, len(chunks), save_frequency)]
                    
                    for batch_idx, chunk_batch in enumerate(chunk_batches):
                        print(f"\n📦 Processing batch {batch_idx + 1}/{len(chunk_batches)} "
                              f"({len(chunk_batch)} chunks)")
                        
                        # Submit batch of chunks
                        future_to_chunk = {
                            executor.submit(self.process_gene_pairs_chunk, chunk, completed_chunks + i): completed_chunks + i 
                            for i, chunk in enumerate(chunk_batch)
                        }
                        
                        # Collect results from this batch
                        batch_results = []
                        for future in as_completed(future_to_chunk):
                            chunk_id = future_to_chunk[future]
                            try:
                                chunk_results = future.result()
                                batch_results.extend(chunk_results)
                                completed_chunks += 1
                                
                                # Progress update
                                total_completed = len(all_results) + len(batch_results)
                                total_to_calc = len(gene_pairs_to_calc)
                                elapsed = time.time() - start_time
                                
                                if total_completed > 0:
                                    eta = (elapsed / total_completed) * (total_to_calc - total_completed)
                                    print(f"✓ Chunk {chunk_id + 1} complete. "
                                          f"Progress: {total_completed:,}/{total_to_calc:,} pairs "
                                          f"({100*total_completed/total_to_calc:.1f}%) - "
                                          f"ETA: {eta/60:.1f} min")
                            
                            except Exception as e:
                                print(f"❌ Chunk {chunk_id + 1} failed: {e}")
                        
                        # Add batch results to total
                        all_results.extend(batch_results)
                        
                        # Save checkpoint after each batch
                        self.save_checkpoint(completed_chunks, total_chunks, all_results)
        
        except KeyboardInterrupt:
            print("\n⚠️ ANALYSIS INTERRUPTED!")
            print("Progress has been saved. You can resume later by running the script again.")
            self.save_checkpoint(completed_chunks, total_chunks, all_results)
            return None
        
        except Exception as e:
            print(f"\n❌ ANALYSIS FAILED: {e}")
            print("Progress has been saved. You can resume after fixing the issue.")
            self.save_checkpoint(completed_chunks, total_chunks, all_results)
            raise
        
        # Create results DataFrame from new calculations
        if all_results:
            new_results_df = pd.DataFrame(all_results, 
                                        columns=['Gene1', 'Gene2', 'jaccard_max', 
                                                'adamic_max', 'gaussian_max'])
        else:
            # Create empty dataframe with correct columns if no new results
            new_results_df = pd.DataFrame(columns=['Gene1', 'Gene2', 'jaccard_max', 
                                                 'adamic_max', 'gaussian_max'])
        
        # Merge with original data
        print(f"\n🔗 Merging results with your original data...")
        
        # If we have existing similarity columns, we need to handle the merge carefully
        if skip_mask is not None:
            # Update only the rows that were calculated
            for _, row in new_results_df.iterrows():
                mask = (input_df['Gene1'] == row['Gene1']) & (input_df['Gene2'] == row['Gene2'])
                input_df.loc[mask, 'jaccard_max'] = row['jaccard_max']
                input_df.loc[mask, 'adamic_max'] = row['adamic_max'] 
                input_df.loc[mask, 'gaussian_max'] = row['gaussian_max']
            
            enhanced_df = input_df
        else:
            # Normal merge for files without existing similarity columns
            enhanced_df = input_df.merge(new_results_df, on=['Gene1', 'Gene2'], how='left')
        
        # DO NOT fill missing values - keep them as None/NaN to preserve real missing data
        # This is the key change from the original code
        
        # Save final results
        print(f"💾 Saving final results to {output_file}...")
        enhanced_df.to_csv(output_file, sep=';', index=False)
        
        # Clean up checkpoint files
        self.cleanup_checkpoint()
        
        total_time = time.time() - start_time
        print(f"\n🎉 ANALYSIS COMPLETE!")
        print(f"   Total time: {total_time/60:.1f} minutes")
        if all_results:
            print(f"   Newly processed: {len(all_results):,} gene pairs")
        if skip_mask is not None:
            skipped = skip_mask.sum()
            print(f"   Skipped (already done): {skipped:,} gene pairs")
        print(f"   Results saved to: {output_file}")
        
        # Print summary statistics
        self.print_summary_stats_enhanced(enhanced_df)
        
        return enhanced_df
    
    def run_analysis(self, output_file: str = "gene_similarity_results.csv", 
                    chunk_size: int = 10000, max_workers: int = None):
        """Run the complete similarity analysis with multi-threading (for all gene pairs)"""
        
        # Load data
        self.load_data()
        
        # Generate all gene pairs
        all_gene_pairs = self.get_all_gene_pairs()
        
        if not all_gene_pairs:
            print("No gene pairs to process!")
            return
        
        # Set default number of workers
        if max_workers is None:
            max_workers = min(os.cpu_count(), 8)  # Don't use too many threads
        
        print(f"Starting analysis with {max_workers} threads, chunk size: {chunk_size:,}")
        
        # Split gene pairs into chunks
        chunks = [all_gene_pairs[i:i + chunk_size] 
                 for i in range(0, len(all_gene_pairs), chunk_size)]
        
        print(f"Created {len(chunks)} chunks")
        
        # Process chunks in parallel
        all_results = []
        start_time = time.time()
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all chunks
            future_to_chunk = {
                executor.submit(self.process_gene_chunk, chunk, i): i 
                for i, chunk in enumerate(chunks)
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_chunk):
                chunk_id = future_to_chunk[future]
                try:
                    chunk_results = future.result()
                    all_results.extend(chunk_results)
                    
                    elapsed = time.time() - start_time
                    completed_pairs = len(all_results)
                    total_pairs = len(all_gene_pairs)
                    
                    if completed_pairs > 0:
                        eta = (elapsed / completed_pairs) * (total_pairs - completed_pairs)
                        print(f"Progress: {completed_pairs:,}/{total_pairs:,} pairs "
                              f"({100*completed_pairs/total_pairs:.1f}%) - "
                              f"ETA: {eta/60:.1f} minutes")
                    
                except Exception as e:
                    print(f"Chunk {chunk_id} failed with error: {e}")
        
        # Save results
        print(f"Saving {len(all_results):,} results to {output_file}...")
        
        results_df = pd.DataFrame(all_results, 
                                columns=['gene_A', 'gene_B', 'jaccard_max', 
                                        'adamic_max', 'gaussian_max'])
        
        results_df.to_csv(output_file, index=False)
        
        total_time = time.time() - start_time
        print(f"Analysis complete! Total time: {total_time/60:.1f} minutes")
        print(f"Results saved to: {output_file}")
        
        # Print summary statistics
        self.print_summary_stats(results_df)
        
    def print_summary_stats_enhanced(self, results_df: pd.DataFrame):
        """Print summary statistics for enhanced results with original columns"""
        print("\n=== ENHANCED SUMMARY STATISTICS ===")
        print(f"Total gene pairs analyzed: {len(results_df):,}")
        
        print(f"\nOriginal columns preserved:")
        original_cols = [col for col in results_df.columns if col not in ['jaccard_max', 'adamic_max', 'gaussian_max']]
        for col in original_cols:
            print(f"  - {col}")
        
        print(f"\nNew similarity methods added:")
        for method in ['jaccard_max', 'adamic_max', 'gaussian_max']:
            if method in results_df.columns:
                values = results_df[method].dropna()  # Only calculate stats for non-null values
                null_count = results_df[method].isna().sum()
                print(f"\n{method.upper()}:")
                if len(values) > 0:
                    print(f"  Mean: {values.mean():.4f}")
                    print(f"  Std:  {values.std():.4f}")
                    print(f"  Min:  {values.min():.4f}")
                    print(f"  Max:  {values.max():.4f}")
                    print(f"  Valid pairs: {len(values):,} ({100*len(values)/len(results_df):.1f}%)")
                    print(f"  Non-zero pairs: {(values > 0).sum():,} ({100*(values > 0).mean():.1f}% of valid)")
                else:
                    print(f"  No valid values calculated")
                print(f"  Missing/None values: {null_count:,} ({100*null_count/len(results_df):.1f}%)")

    def print_summary_stats(self, results_df: pd.DataFrame):
        """Print summary statistics for standard results"""
        print("\n=== SUMMARY STATISTICS ===")
        print(f"Total gene pairs analyzed: {len(results_df):,}")
        
        for method in ['jaccard_max', 'adamic_max', 'gaussian_max']:
            if method in results_df.columns:
                values = results_df[method].dropna()  # Only calculate stats for non-null values
                null_count = results_df[method].isna().sum()
                print(f"\n{method.upper()}:")
                if len(values) > 0:
                    print(f"  Mean: {values.mean():.4f}")
                    print(f"  Std:  {values.std():.4f}")
                    print(f"  Min:  {values.min():.4f}")
                    print(f"  Max:  {values.max():.4f}")
                    print(f"  Valid pairs: {len(values):,} ({100*len(values)/len(results_df):.1f}%)")
                    print(f"  Non-zero pairs: {(values > 0).sum():,} ({100*(values > 0).mean():.1f}% of valid)")
                else:
                    print(f"  No valid values calculated")
                print(f"  Missing/None values: {null_count:,} ({100*null_count/len(results_df):.1f}%)")

def main():
    """Main execution function - Process existing gene pairs with checkpoint support"""
    
    # Configuration
    DATA_DIR = "./"  # Directory containing depth_X.csv files
    INPUT_CSV = "HSFinalSIVF.csv"  # Your existing gene pairs file
    OUTPUT_FILE = "HSFinalSIVF_enhanced.csv"  # Output with added columns
    CHUNK_SIZE = 2000  # Smaller chunks for more frequent checkpoints
    MAX_WORKERS = 6    # Adjust based on your CPU cores
    SAVE_FREQUENCY = 3  # Save checkpoint every 3 chunks
    SIGMA = 1.2        # Gaussian kernel parameter (optimized for binary vectors)
    MIN_FREQUENCY = 2  # Minimum frequency for Adamic/Adar (avoids log(1)=0)
    
    print("🧬 Gene Similarity Analyzer with Checkpoint Support (Keep None Values)")
    print("=" * 70)
    
    # Create analyzer
    analyzer = GeneSimilarityAnalyzer(
        data_dir=DATA_DIR,
        sigma=SIGMA,
        min_frequency=MIN_FREQUENCY
    )
    
    # Debug data structure first
    print("=== DEBUGGING DATA STRUCTURE ===")
    analyzer.debug_data_structure()
    
    # Process your existing gene pairs file
    if os.path.exists(INPUT_CSV):
        print(f"📁 Found your gene pairs file: {INPUT_CSV}")
        print("🔬 Adding similarity calculations to your existing data...")
        print(f"⚙️ Configuration:")
        print(f"   - Chunk size: {CHUNK_SIZE:,} pairs")
        print(f"   - Threads: {MAX_WORKERS}")
        print(f"   - Auto-save: Every {SAVE_FREQUENCY} chunks")
        print(f"   - Sigma: {SIGMA}")
        print(f"   - Min frequency: {MIN_FREQUENCY}")
        print(f"   - Missing genes: Will use None (not 0.0)")
        
        try:
            enhanced_df = analyzer.run_analysis_from_csv(
                input_csv_file=INPUT_CSV,
                output_file=OUTPUT_FILE,
                chunk_size=CHUNK_SIZE,
                max_workers=MAX_WORKERS,
                save_frequency=SAVE_FREQUENCY
            )
            
            if enhanced_df is not None:
                print(f"\n📊 Sample of enhanced results:")
                print(enhanced_df.head())
            
        except Exception as e:
            print(f"💥 Analysis failed with error: {e}")
            print("🔄 Progress has been saved - you can resume by running the script again.")
            raise
    else:
        print(f"❌ ERROR: Your gene pairs file '{INPUT_CSV}' not found!")
        print(f"📍 Please make sure '{INPUT_CSV}' is in the same directory as this script.")
        print("📋 Expected format:")
        print("   Gene1;Gene2;Interaction;Similarity")
        print("   UBE2Q1;RNF114;902;0.6792731330979971")
        print("   UBE2Q1;FBXL12;900;0.5458068064294364")

if __name__ == "__main__":
    main()