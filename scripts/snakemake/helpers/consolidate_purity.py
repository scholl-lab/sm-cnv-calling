#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import sys

def get_final_purity(purecn_csv, fallback_purity, output_file):
    try:
        assert 0.0 < fallback_purity <= 1.0, f"Fallback purity must be between 0 and 1. Got: {fallback_purity}"
    except AssertionError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
        
    final_purity = fallback_purity
    source = "fallback"

    if purecn_csv and os.path.exists(purecn_csv):
        try:
            df = pd.read_csv(purecn_csv)
            if not df.empty and 'Purity' in df.columns:
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
    parser = argparse.ArgumentParser(description="Consolidate tumor purity estimates.")
    parser.add_argument('--purecn-csv', required=False, help="Path to the PureCN output CSV file. Can be 'None'.")
    parser.add_argument('--fallback-purity', type=float, required=True, help="Fallback purity estimate.")
    parser.add_argument('--output', required=True, help="Path to write the final purity value.")
    args = parser.parse_args()
    purecn_path = None if args.purecn_csv == 'None' else args.purecn_csv
    get_final_purity(purecn_path, args.fallback_purity, args.output)
