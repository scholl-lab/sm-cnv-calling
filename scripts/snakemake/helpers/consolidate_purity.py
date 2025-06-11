import argparse
import pandas as pd
import os

def get_final_purity(purecn_csv, fallback_purity, output_file):
    """
    Determines the final purity value.
    Prioritizes the top hit from a PureCN CSV file if it exists and is valid.
    Otherwise, uses the provided fallback value.
    """
    final_purity = fallback_purity
    source = "fallback"

    if purecn_csv and os.path.exists(purecn_csv):
        try:
            df = pd.read_csv(purecn_csv)
            if not df.empty and 'Purity' in df.columns:
                # Get the purity from the first row (top solution)
                top_purity = df.iloc[0]['Purity']
                if pd.notna(top_purity) and 0 < top_purity <= 1.0:
                    final_purity = top_purity
                    source = "PureCN"
                else:
                    print(f"Warning: Invalid purity value '{top_purity}' in {purecn_csv}. Using fallback.")
            else:
                print(f"Warning: PureCN output {purecn_csv} is empty or missing 'Purity' column. Using fallback.")
        except Exception as e:
            print(f"Warning: Could not parse {purecn_csv}. Error: {e}. Using fallback.")
    
    print(f"Final purity for {os.path.basename(output_file).replace('.purity.txt','')} is {final_purity} (source: {source})")
    
    with open(output_file, 'w') as f_out:
        f_out.write(str(final_purity))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Consolidate tumor purity estimates, using PureCN output with a fallback.")
    parser.add_argument('--purecn-csv', required=False, help="Path to the PureCN output CSV file. Can be 'None' if not applicable.")
    parser.add_argument('--fallback-purity', type=float, required=True, help="Fallback purity estimate if PureCN is not available or fails.")
    parser.add_argument('--output', required=True, help="Path to write the final purity value.")
    
    args = parser.parse_args()

    # Handle the 'None' string from Snakemake params
    purecn_path = None if args.purecn_csv == 'None' else args.purecn_csv

    get_final_purity(purecn_path, args.fallback_purity, args.output)
