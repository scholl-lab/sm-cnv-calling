#!/usr/bin/env python3
import argparse
import os
import sys
import yaml
import requests
import gzip
from tqdm import tqdm

def download_file(url, output_path):
    """Downloads a file from a URL to a given path with a progress bar."""
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for bad status codes
        total_size = int(response.headers.get('content-length', 0))
        
        # Check if the URL points to a gzipped file
        is_gzipped = url.endswith('.gz')
        final_output_path = output_path.replace('.gz', '') if is_gzipped else output_path

        print(f"Downloading {os.path.basename(url)} to {final_output_path}...")
        
        with open(final_output_path, 'wb') as f, tqdm(
            desc=os.path.basename(final_output_path),
            total=total_size,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            # If gzipped, decompress on the fly
            if is_gzipped:
                # requests.iter_content doesn't work well with gzip decompression on the fly,
                # so we download to a temporary file first.
                temp_gzipped_path = final_output_path + ".gz.tmp"
                with open(temp_gzipped_path, 'wb') as f_tmp:
                    for chunk in response.iter_content(chunk_size=8192):
                        bar.update(len(chunk))
                        f_tmp.write(chunk)
                
                # Now decompress from the temp file
                with gzip.open(temp_gzipped_path, 'rb') as f_gz:
                    f.write(f_gz.read())
                os.remove(temp_gzipped_path)
            else:
                 for chunk in response.iter_content(chunk_size=8192):
                    bar.update(len(chunk))
                    f.write(chunk)
        print(f"Successfully downloaded and saved {final_output_path}\n")

    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Download standard accessory files for a reference genome.")
    parser.add_argument('-g', '--genome', required=True, choices=['hg19', 'hg38'], help="Reference genome build (e.g., 'hg19', 'hg38').")
    parser.add_argument('-o', '--output-dir', required=True, help="Directory to save the downloaded files.")
    parser.add_argument('-c', '--config', default='download_config.yaml', help="Path to the downloader config YAML file.")
    
    args = parser.parse_args()

    # Load configuration
    try:
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Config file not found at {args.config}", file=sys.stderr)
        sys.exit(1)

    if args.genome not in config:
        print(f"Error: Genome build '{args.genome}' not found in config file {args.config}", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    genome_urls = config[args.genome]
    
    # Download access bed file
    if 'access_bed' in genome_urls:
        url = genome_urls['access_bed']
        filename = os.path.basename(url)
        download_file(url, os.path.join(args.output_dir, filename))
        
    # Download and decompress refFlat file
    if 'ref_flat' in genome_urls:
        url = genome_urls['ref_flat']
        # The output filename will be refFlat.txt (without .gz)
        filename = os.path.basename(url)
        download_file(url, os.path.join(args.output_dir, filename))

    print("All downloads complete.")

if __name__ == '__main__':
    main()
