import pandas as pd
import os

def analyze_csv_files():
    # Get list of all CSV files in current directory
    csv_files = [f for f in os.listdir('.') if f.endswith('.csv')]
    
    results = {}
    
    for file in csv_files:
        try:
            # Read CSV file
            df = pd.read_csv(file)
            
            # Get number of rows and columns
            num_rows = len(df)
            num_cols = len(df.columns)
            
            # Store results
            results[file] = {
                'rows': num_rows,
                'columns': num_cols,
                'shape': df.shape
            }
            
            print(f"\nAnalyzing {file}:")
            print(f"Number of rows: {num_rows}")
            print(f"Number of columns: {num_cols}")
            print(f"Shape: {df.shape}")
            
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
    
    return results

if __name__ == "__main__":
    results = analyze_csv_files()
