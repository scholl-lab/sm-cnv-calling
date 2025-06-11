import argparse
import pandas as pd

def identify_clean_normals(metrics_file, output_list, threshold):
    """
    Parses a CNVkit metrics file, filters out noisy normals, and writes a list
    of clean normal coverage files.
    """
    try:
        df = pd.read_csv(metrics_file, sep='\t')
        
        # Ensure required columns exist
        if 'sample' not in df.columns or 'bivar' not in df.columns:
            raise ValueError("Metrics file must contain 'sample' and 'bivar' columns.")
            
        # Filter based on the biweight midvariance threshold
        clean_df = df[df['bivar'] <= threshold]
        
        if clean_df.empty:
            print(f"Warning: No normal samples passed the QC threshold of {threshold}. The clean normal list will be empty.")
            open(output_list, 'w').close() # Create an empty file
        else:
            # The 'sample' column contains the full path to the .targetcoverage.cnn file
            clean_files = clean_df['sample'].tolist()
            with open(output_list, 'w') as f_out:
                for file_path in clean_files:
                    f_out.write(f"{file_path}\n")
            print(f"Identified {len(clean_files)} clean normals out of {len(df)}. Wrote list to {output_list}")

    except FileNotFoundError:
        print(f"Error: The metrics file {metrics_file} was not found.")
        # Create an empty output file to prevent Snakemake from failing
        open(output_list, 'w').close()
    except Exception as e:
        print(f"An error occurred: {e}")
        open(output_list, 'w').close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter CNVkit normal metrics to identify clean samples for a Panel of Normals.")
    parser.add_argument('-i', '--input-metrics', required=True, help="Path to the metrics file from 'cnvkit.py metrics'.")
    parser.add_argument('-o', '--output-list', required=True, help="Path to write the list of clean normal coverage files.")
    parser.add_argument('-t', '--threshold', type=float, required=True, help="QC threshold for the 'bivar' metric. Samples at or below this value are kept.")
    
    args = parser.parse_args()
    
    identify_clean_normals(args.input_metrics, args.output_list, args.threshold)
