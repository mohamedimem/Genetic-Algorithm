from step3_commondepath import complete_gene_analysis_standalone
import csv
import time

def analyze_gene_pairs():
    results = []
    
    # First, read all rows into a list
    all_rows = []
    with open('HSFinalSIVF.csv', mode='r') as infile:
        reader = csv.reader(infile, delimiter=';')
        next(reader)  # Skip header row
        all_rows = list(reader)  # Read all remaining rows into a list
    
    # Now process the rows in reverse order (last to first)
    for row in reversed(all_rows):
        gene1 = row[0]
        gene2 = row[1]
        try:
            location_score = complete_gene_analysis_standalone(gene1, gene2)
            results.append((gene1, gene2, location_score))
            print(f"Processed: {gene1}-{gene2}, Score: {location_score}")  # Optional: track progress
          #  time.sleep(1)  # Add 1 second delay between analyses

        except Exception as e:
            print(f"Error analyzing {gene1}-{gene2}: {e}")
            continue
    
    return results

# Call the function
results = analyze_gene_pairs()
print(f"Analysis complete. Processed {len(results)} gene pairs.")