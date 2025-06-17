# In your other module
from step3_commondepath import (
    GeneCommonAncestorAnalyzer, 
    normalize_depth_score,
    complete_gene_analysis_standalone  # if you want the complete analysis too
)

def analyze_gene_similarity(gene1, gene2):
    """Your custom analysis function"""
    analyzer = GeneCommonAncestorAnalyzer()
    score = analyzer.complete_gene_analysis(gene1, gene2, verbose=False)
    
    return score

# Usage

test_pairs = [
        ("CD19", "ACTG1"),      # Actin cytoskeleton components
        ("TUBB", "TUBA1A"),     # Tubulin/microtubule components  
        ("COX1", "COX2"),       # Mitochondrial respiratory chain complex IV
        ("ATP5F1A", "ATP5F1B"), # ATP synthase complex components
        ("HIST1H1A", "HIST1H2A"), # Histone/chromatin components
    ]

    

for gene1, gene2 in test_pairs:
    result = analyze_gene_similarity(gene1, gene2)
    print(f" {gene1} vs {gene2} Location score: {result:.2f}")