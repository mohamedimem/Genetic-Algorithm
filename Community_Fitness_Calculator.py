import csv
import sys
from datetime import datetime

def create_Dictionnairefinal(): 
    """Load gene interaction and similarity data from CSV"""
    l1 = list()
    l2 = list()
    
    try:
        with open('HSFinalSIVF.csv', mode='r') as infile:
            reader = csv.reader(infile, delimiter=';')
            next(reader)  # Skip header
            for rows in reader:
                g = (rows[0], rows[1])                
                l1.append(g) 
                IS = (float(rows[2])/1000, float(rows[3]))  
                l2.append(IS) 
                
        DictF = {k:v for k, v in zip(l1, l2)}
        print(f"Loaded {len(DictF)} gene pairs from HSFinalSIVF.csv")
        return DictF
        
    except FileNotFoundError:
        print("ERROR: HSFinalSIVF.csv file not found!")
        print("Please ensure the CSV file is in the same directory as this script.")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR reading CSV file: {e}")
        sys.exit(1)

def simInt_fsi(community, dictgene, verbose=True): 
    """
    Calculate fitness for a gene community using the original GA fitness function
    Returns: (avg_interaction, avg_similarity, fitness_score)
    """
    st = 0  # Total number of pairs
    si = 0.0  # Sum of interactions
    ss = 0.0  # Sum of similarities
    alpha = 0.5 
    beta = 0.5 
    
    pairs_found = 0
    pairs_missing = 0
    pair_details = []
    
    if verbose:
        print(f"\nAnalyzing community of {len(community)} genes:")
        print(f"Genes: {', '.join(community)}")
        print(f"\nPair-by-pair analysis:")
        print("-" * 80)
    
    for i, item in enumerate(community):
        for j in range(i+1, len(community)):
            g1 = community[i]
            g2 = community[j]
            genes = (g1, g2)
            geness = (g2, g1)
            t = len(community)
            
            v1 = dictgene.get(genes, 0) 
            v2 = dictgene.get(geness, 0) 
            
            if v1 == 0: 
                v = v2 
            else: 
                v = v1 
                
            if v != 0:
                intr = v[0]
                s = v[1]
                si = si + intr 
                ss = ss + s
                pairs_found += 1
                
                if verbose:
                    print(f"{g1} - {g2}: Interaction={intr:.6f}, Similarity={s:.6f} [FOUND]")
                
                pair_details.append({
                    'pair': f"{g1}-{g2}",
                    'interaction': intr,
                    'similarity': s,
                    'found': True
                })
                
            else:
                intr = 0.0
                s = 0.0
                pairs_missing += 1
                
                if verbose:
                    print(f"{g1} - {g2}: No data found (using 0.0, 0.0) [MISSING]")
                
                pair_details.append({
                    'pair': f"{g1}-{g2}",
                    'interaction': 0.0,
                    'similarity': 0.0,
                    'found': False
                })
                
        st = st + (t-i)
    
    # Calculate averages
    sint = si/st
    smoyen = ss/st
    
    # Calculate fitness using original formula
    f = alpha*(smoyen) + beta*(sint)
    
    if verbose:
        print("-" * 80)
        print(f"SUMMARY STATISTICS:")
        print(f"Total gene pairs: {st}")
        print(f"Pairs found in CSV: {pairs_found}")
        print(f"Pairs missing from CSV: {pairs_missing}")
        print(f"Coverage: {(pairs_found/st)*100:.1f}%")
        print(f"\nCOMPONENT SCORES:")
        print(f"Average Interaction: {sint:.6f}")
        print(f"Average Similarity: {smoyen:.6f}")
        print(f"\nFITNESS CALCULATION:")
        print(f"Formula: {alpha} × AVG_SIM + {beta} × AVG_INT")
        print(f"Calculation: {alpha} × {smoyen:.6f} + {beta} × {sint:.6f}")
        print(f"FINAL FITNESS: {f:.6f}")
    
    return sint, smoyen, f, pair_details

def analyze_gene_contributions(community, dictgene):
    """Analyze individual gene contributions to the community fitness"""
    print(f"\nINDIVIDUAL GENE ANALYSIS:")
    print("-" * 60)
    
    gene_scores = []
    
    for index, target_gene in enumerate(community):
        si = 0.0
        ss = 0.0
        connections = 0
        
        for j, other_gene in enumerate(community):
            if index != j:  # Don't compare gene with itself
                genes = (target_gene, other_gene)
                geness = (other_gene, target_gene)
                
                v1 = dictgene.get(genes, 0) 
                v2 = dictgene.get(geness, 0) 
                
                if v1 == 0: 
                    v = v2 
                else: 
                    v = v1 
                    
                if v != 0:
                    i = float(v[0])    
                    s = float(v[1])
                    connections += 1
                else: 
                    i = 0.0
                    s = 0.0
                    
                si += i
                ss += s
        
        # Calculate average scores for this gene
        t = len(community) - 1  # Exclude self
        avg_sim = ss / t if t > 0 else 0
        avg_int = si / t if t > 0 else 0
        gene_score = avg_sim + avg_int
        
        gene_scores.append({
            'gene': target_gene,
            'avg_similarity': avg_sim,
            'avg_interaction': avg_int,
            'score': gene_score,
            'connections': connections
        })
        
        print(f"{target_gene:8s}: Score={gene_score:.6f} (Sim={avg_sim:.6f}, Int={avg_int:.6f}, Connections={connections}/{t})")
    
    # Sort by score
    gene_scores.sort(key=lambda x: x['score'], reverse=True)
    
    print(f"\nGENE RANKING BY CONTRIBUTION:")
    print("-" * 60)
    for i, gene_data in enumerate(gene_scores, 1):
        print(f"{i:2d}. {gene_data['gene']:8s} - Score: {gene_data['score']:.6f}")
    
    return gene_scores

def main():
    """Main function to calculate fitness for the given community"""
    
    print("=" * 80)
    print("GENETIC ALGORITHM FITNESS CALCULATOR")
    print("=" * 80)
    print(f"Analysis started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Define the community to analyze
     
    given_community = ['TXNL4A', 'DDX39A', 'SRSF9', 'UPF3B', 'SSRP1', 'POLR2I', 'SRSF6']
    
    print(f"\nTarget Community: {', '.join(given_community)}")
    print(f"Community Size: {len(given_community)} genes")
    
    # Load the gene interaction data
    dictgene = create_Dictionnairefinal()
    
    # Calculate fitness
    sint, smoyen, fitness, pair_details = simInt_fsi(given_community, dictgene, verbose=True)
    
    # Analyze individual gene contributions
    gene_scores = analyze_gene_contributions(given_community, dictgene)
    
    # Save detailed results to file
    output_file = f"Community_Fitness_Analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("GENETIC ALGORITHM FITNESS ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Community: {', '.join(given_community)}\n")
        f.write(f"Size: {len(given_community)} genes\n\n")
        
        f.write("FITNESS RESULTS:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Average Interaction: {sint:.6f}\n")
        f.write(f"Average Similarity:  {smoyen:.6f}\n")
        f.write(f"Final Fitness:       {fitness:.6f}\n\n")
        
        f.write("PAIR-BY-PAIR BREAKDOWN:\n")
        f.write("-" * 40 + "\n")
        for detail in pair_details:
            status = "[FOUND]" if detail['found'] else "[MISSING]"
            f.write(f"{detail['pair']:15s}: Int={detail['interaction']:.6f}, Sim={detail['similarity']:.6f} {status}\n")
        
        f.write(f"\nGENE RANKING:\n")
        f.write("-" * 40 + "\n")
        for i, gene_data in enumerate(gene_scores, 1):
            f.write(f"{i:2d}. {gene_data['gene']:8s} - Score: {gene_data['score']:.6f}\n")
    
    print(f"\n" + "=" * 80)
    print(f"FINAL RESULTS:")
    print(f"Community Fitness: {fitness:.6f}")
    print(f"Best Contributing Gene: {gene_scores[0]['gene']} (Score: {gene_scores[0]['score']:.6f})")
    print(f"Detailed results saved to: {output_file}")
    print("=" * 80)

if __name__ == '__main__':
    main()