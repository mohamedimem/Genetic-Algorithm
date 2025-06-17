import pandas as pd
import numpy as np
import os
from typing import List, Tuple, Optional, Dict
from uniport_step1 import UniProtGOExtractor
from step2_geolocation import download_go_obo, GODag, get_go_info, save_gene_depths_to_csv

# Global shared instances to prevent repeated loading
_SHARED_GO_DAG = None
_SHARED_EXTRACTOR = None
_SHARED_COUNT_DF = None
_COUNT_DF_LAST_MODIFIED = None

def get_shared_go_dag():
    """Get or create the shared GO DAG instance."""
    global _SHARED_GO_DAG
    if _SHARED_GO_DAG is None:
        print("🔄 Loading GO DAG (one-time initialization)...")
        _SHARED_GO_DAG = GODag("go-basic.obo")
        print("✅ GO DAG loaded and cached")
    return _SHARED_GO_DAG

def get_shared_extractor():
    """Get or create the shared UniProt extractor instance."""
    global _SHARED_EXTRACTOR
    if _SHARED_EXTRACTOR is None:
        print("🔄 Initializing UniProt extractor (one-time initialization)...")
        _SHARED_EXTRACTOR = UniProtGOExtractor()
        print("✅ UniProt extractor initialized and cached")
    return _SHARED_EXTRACTOR

def get_shared_count_df(count_depth_file="count_depth.csv", force_reload=False):
    """Get or create the shared count depth DataFrame with smart caching."""
    global _SHARED_COUNT_DF, _COUNT_DF_LAST_MODIFIED
    
    if not os.path.exists(count_depth_file):
        raise FileNotFoundError(f"Count depth file '{count_depth_file}' not found!")
    
    # Check if file has been modified since last load
    current_mtime = os.path.getmtime(count_depth_file)
    
    if force_reload or _SHARED_COUNT_DF is None or _COUNT_DF_LAST_MODIFIED != current_mtime:
        print(f"🔄 Reloading count depth data (file updated): {len(pd.read_csv(count_depth_file, index_col=0))} genes")
        _SHARED_COUNT_DF = pd.read_csv(count_depth_file, index_col=0)
        _COUNT_DF_LAST_MODIFIED = current_mtime
    
    return _SHARED_COUNT_DF

class GeneCommonAncestorAnalyzer:
    """
    A class to analyze common ancestors between genes using depth-based CSV files
    and Jaccard similarity on binary vectors.
    Uses global shared instances to prevent repeated loading.
    """
    
    def __init__(self, count_depth_file="count_depth.csv"):
        """
        Initialize the analyzer.
        
        Args:
            count_depth_file (str): Path to the count_depth.csv file
        """
        self.count_depth_file = count_depth_file
        # Use shared instances instead of creating new ones
        self.count_df = get_shared_count_df(count_depth_file)
        self.go_dag = get_shared_go_dag()
        self.extractor = get_shared_extractor()
    
    def _reload_count_depth_if_needed(self):
        """Reload count depth data only if the file has been modified."""
        self.count_df = get_shared_count_df(self.count_depth_file, force_reload=False)
    
    def gene_exists_in_count_depth(self, gene_name: str) -> bool:
        """
        Check if a gene exists in the count_depth.csv file.
        
        Args:
            gene_name (str): Name of the gene to check
            
        Returns:
            bool: True if gene exists in count_depth data, False otherwise
        """
        if self.count_df is None:
            print("Warning: Count depth data not loaded")
            return False
        
        return gene_name in self.count_df.index
    
    def get_common_depth_levels(self, gene1: str, gene2: str) -> List[int]:
        """
        Find depth levels where both genes have non-zero counts.
        
        Args:
            gene1 (str): First gene name
            gene2 (str): Second gene name
            
        Returns:
            List[int]: List of common depth levels, sorted from highest to lowest
        """
        if gene1 not in self.count_df.index:
            return []
        
        if gene2 not in self.count_df.index:
            return []
        
        gene1_data = self.count_df.loc[gene1]
        gene2_data = self.count_df.loc[gene2]
        
        common_depths = []
        
        # Check each depth column
        for col in self.count_df.columns:
            if col.startswith('depth_') and col.endswith('_count'):
                # Extract depth number from column name (e.g., 'depth_3_count' -> 3)
                try:
                    depth_num = int(col.split('_')[1])
                    
                    # Check if both genes have non-zero counts at this depth
                    gene1_count = gene1_data[col] if not pd.isna(gene1_data[col]) else 0
                    gene2_count = gene2_data[col] if not pd.isna(gene2_data[col]) else 0
                    
                    if gene1_count > 0 and gene2_count > 0:
                        common_depths.append(depth_num)
                        
                except (ValueError, IndexError):
                    continue
        
        # Sort from highest to lowest (11 -> 10 -> 9 -> ...)
        return sorted(common_depths, reverse=True)
    
    def load_depth_csv(self, depth_level: int) -> Optional[pd.DataFrame]:
        """
        Load a specific depth CSV file.
        
        Args:
            depth_level (int): The depth level to load
            
        Returns:
            pd.DataFrame or None: The loaded DataFrame or None if file doesn't exist
        """
        depth_file = f"depth_{depth_level}.csv"
        
        if not os.path.exists(depth_file):
            return None
        
        try:
            df = pd.read_csv(depth_file, index_col=0)
            return df
        except Exception as e:
            print(f"Error loading depth file '{depth_file}': {e}")
            return None
    
    def get_gene_binary_vectors(self, depth_df: pd.DataFrame, gene1: str, gene2: str) -> Tuple[Optional[np.array], Optional[np.array]]:
        """
        Extract binary vectors for two genes from a depth DataFrame.
        
        Args:
            depth_df (pd.DataFrame): The depth DataFrame
            gene1 (str): First gene name
            gene2 (str): Second gene name
            
        Returns:
            Tuple[np.array, np.array]: Binary vectors for both genes (or None if gene not found)
        """
        vector1 = None
        vector2 = None
        
        if gene1 in depth_df.index:
            vector1 = depth_df.loc[gene1].values.astype(int)
        
        if gene2 in depth_df.index:
            vector2 = depth_df.loc[gene2].values.astype(int)
        
        return vector1, vector2
    
    def calculate_jaccard_similarity(self, vector1: np.array, vector2: np.array) -> float:
        """
        Calculate Jaccard similarity between two binary vectors.
        
        Args:
            vector1 (np.array): First binary vector
            vector2 (np.array): Second binary vector
            
        Returns:
            float: Jaccard similarity coefficient (0 to 1)
        """
        if vector1 is None or vector2 is None:
            return 0.0
        
        if len(vector1) != len(vector2):
            return 0.0
        
        # Convert to boolean arrays for intersection and union
        v1_bool = vector1.astype(bool)
        v2_bool = vector2.astype(bool)
        
        # Calculate intersection and union
        intersection = np.sum(v1_bool & v2_bool)
        union = np.sum(v1_bool | v2_bool)
        
        # Avoid division by zero
        if union == 0:
            return 0.0
        
        jaccard = intersection / union
        return jaccard
    
    def find_deepest_nonzero_jaccard(self, gene1: str, gene2: str, verbose: bool = False) -> Tuple[float, Optional[int]]:
        """
        Find Jaccard similarity from the DEEPEST level that has NON-ZERO similarity.
        
        Returns:
            Tuple[float, Optional[int]]: (Jaccard similarity, depth level) from deepest non-zero level
        """
        common_depths = self.get_common_depth_levels(gene1, gene2)
        if not common_depths:
            return 0.0, None
        
        # Process from deepest to shallowest
        for depth in common_depths:
            depth_df = self.load_depth_csv(depth)
            if depth_df is None:
                continue
            
            vector1, vector2 = self.get_gene_binary_vectors(depth_df, gene1, gene2)
            if vector1 is None or vector2 is None:
                continue
            
            jaccard_score = self.calculate_jaccard_similarity(vector1, vector2)
            
            # Return first NON-ZERO score (which is the deepest)
            if jaccard_score > 0:
                return jaccard_score, depth
        
        return 0.0, None
    
    def complete_gene_analysis(self, gene1: str, gene2: str, verbose: bool = False) -> float:
        """
        Complete analysis workflow for two genes including:
        1. Check if genes exist in count_depth.csv
        2. If not, fetch GO annotations and create depth representations
        3. Analyze gene pair similarity
        
        Args:
            gene1 (str): First gene name
            gene2 (str): Second gene name
            verbose (bool): Whether to print detailed progress information
            
        Returns:
            float: Location score (0.0 if error)
        """
        # Check if genes exist in count_depth.csv
        gene1_exists = self.gene_exists_in_count_depth(gene1)
        gene2_exists = self.gene_exists_in_count_depth(gene2)
        
        # Process genes that don't exist in count_depth.csv
        genes_to_process = []
        if not gene1_exists:
            genes_to_process.append(gene1)
        if not gene2_exists:
            genes_to_process.append(gene2)
        
        if genes_to_process:
            if verbose:
                print(f"📝 Processing missing genes: {genes_to_process}")
            
            for gene in genes_to_process:
                try:
                    # Step 1: Fetch GO annotations (using shared extractor)
                    go_cc_terms = self.extractor.step1_fetch_go_annotations(gene)
                    
                    if not go_cc_terms:
                        if verbose:
                            print(f"⚠️ No GO annotations found for {gene}")
                        continue
                    
                    # Step 2: Get GO info and create depth representations (using shared GO DAG)
                    target_go_ids = list(go_cc_terms)
                    go_info_results = get_go_info(target_go_ids, self.go_dag)
                    
                    # Step 3: Save gene depths to CSV
                    save_gene_depths_to_csv(gene, go_info_results)
                    if verbose:
                        print(f"✅ Processed {gene} successfully")
                    
                except Exception as e:
                    if verbose:
                        print(f"❌ Error processing {gene}: {e}")
                    return 0.0
            
            # IMPORTANT: Reload count_depth data after processing new genes
            # This ensures we can see the newly added gene vectors
            if verbose:
                print("🔄 Reloading count depth data after gene processing...")
            self._reload_count_depth_if_needed()
        
        # Final check that both genes now exist
        gene1_exists_final = self.gene_exists_in_count_depth(gene1)
        gene2_exists_final = self.gene_exists_in_count_depth(gene2)
        
        if not gene1_exists_final or not gene2_exists_final:
            missing_genes = []
            if not gene1_exists_final:
                missing_genes.append(gene1)
            if not gene2_exists_final:
                missing_genes.append(gene2)
            
            if verbose:
                print(f"❌ Genes still missing after processing: {missing_genes}")
            return 0.0
        
        # Step 4: Analyze gene pair
        try:
            #score, depth = self.find_deepest_nonzero_jaccard(gene1, gene2, verbose=False)
           # location_score = normalize_depth_score(depth)
            location_score=0
            return location_score
            
        except Exception as e:
            if verbose:
                print(f"❌ Error during gene pair analysis: {e}")
            return 0.0

# Helper functions
def normalize_depth_score(depth):
    """
    Normalize depth value to location score.
    
    Args:
        depth (int): Depth level
        
    Returns:
        float: Normalized score (1.0, 0.8, or 0.6)
    """
    if depth is None:
        return 0.0
    elif depth > 5:
        return 1.0
    elif depth > 3:
        return 0.8
    else:
        return 0.6

def reset_shared_cache():
    """Reset all shared cached objects (for testing or cleanup)."""
    global _SHARED_GO_DAG, _SHARED_EXTRACTOR, _SHARED_COUNT_DF, _COUNT_DF_LAST_MODIFIED
    _SHARED_GO_DAG = None
    _SHARED_EXTRACTOR = None
    _SHARED_COUNT_DF = None
    _COUNT_DF_LAST_MODIFIED = None
    print("🧹 Shared cache reset")

# Simplified usage functions
def complete_gene_analysis_standalone(gene1: str, gene2: str, verbose: bool = False) -> float:
    """
    Standalone version of complete gene analysis that uses shared instances.
    
    Args:
        gene1 (str): First gene name
        gene2 (str): Second gene name
        verbose (bool): Whether to print detailed progress information
        
    Returns:
        float: Location score (0.0 if error or no similarity found)
    """
    try:
        # This will reuse shared instances automatically
        analyzer = GeneCommonAncestorAnalyzer()
        result = analyzer.complete_gene_analysis(gene1, gene2, verbose)
        return float(result)
            
    except Exception as e:
        if verbose:
            print(f"❌ Exception in analysis: {e}")
        return 0.0

# Additional utility function for batch processing
def create_persistent_analyzer():
    """
    Create an analyzer instance that will persist and reuse shared resources.
    Use this when you plan to analyze many gene pairs in sequence.
    
    Returns:
        GeneCommonAncestorAnalyzer: Analyzer with persistent shared resources
    """
    print("🚀 Creating persistent analyzer with shared resources...")
    
    # Pre-load shared resources
    get_shared_go_dag()
    get_shared_extractor() 
    get_shared_count_df()
    
    analyzer = GeneCommonAncestorAnalyzer()
    print("✅ Persistent analyzer ready for batch processing")
    return analyzer

def batch_gene_analysis(gene_pairs: List[Tuple[str, str]], verbose: bool = False) -> List[float]:
    """
    Efficiently analyze multiple gene pairs using a single persistent analyzer.
    
    Args:
        gene_pairs: List of (gene1, gene2) tuples to analyze
        verbose: Whether to print progress
        
    Returns:
        List[float]: Location scores for each gene pair
    """
    analyzer = create_persistent_analyzer()
    results = []
    
    for i, (gene1, gene2) in enumerate(gene_pairs):
        if verbose:
            print(f"\n📊 Analyzing pair {i+1}/{len(gene_pairs)}: {gene1} vs {gene2}")
        
        score = analyzer.complete_gene_analysis(gene1, gene2, verbose=False)
        results.append(score)
        
        if verbose:
            print(f"   Result: {score:.2f}")
    
    return results