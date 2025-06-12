import os
import sys
import pandas as pd

# ----------------------------------------------------------------------------------- #
#                         Main Controller for CNV-Pipline                             #
# ----------------------------------------------------------------------------------- #

# --- Configuration and Setup ---
configfile: "config.yaml"

# --- Path Validation ---
# Fail fast if essential reference files are missing
for ref_key in ["reference_genome", "targets_bed", "access_bed"]:
    if not os.path.exists(config[ref_key]):
        sys.exit(f"Critical Error: Reference file not found at path specified by '{ref_key}': {config[ref_key]}")

# Validate optional annotation file if specified
if config.get("annotate_refFlat") and not os.path.exists(config["annotate_refFlat"]):
    sys.exit(f"Critical Error: Optional annotation file not found at path specified by 'annotate_refFlat': {config['annotate_refFlat']}")

# --- Sample Sheet Loading and Helpers ---
try:
    SAMPLES = pd.read_csv(config["samplesheet"], sep="\t", dtype=str)
    # Check for duplicate sample IDs
    if SAMPLES["sample_id"].duplicated().any():
        duplicated_ids = SAMPLES[SAMPLES["sample_id"].duplicated(keep=False)]["sample_id"].tolist()
        sys.exit(f"Critical Error: Duplicate sample IDs found in sample sheet: {duplicated_ids}")
    SAMPLES = SAMPLES.set_index("sample_id", drop=False)
except FileNotFoundError:
    sys.exit(f"Critical Error: Sample sheet not found at {config['samplesheet']}. Please create it.")

def get_tumor_samples():
    return SAMPLES[SAMPLES.analysis_type.isin(["TvsN", "To"])].index.tolist()

def get_normal_bams_for_pon():
    """Get deduplicated normal BAM files from TvsN samples for Panel of Normals construction"""
    tvsn_samples = SAMPLES[SAMPLES.analysis_type == "TvsN"]
    normal_bams = {}
    seen_bams = set()
    
    for sample_id, row in tvsn_samples.iterrows():
        if pd.notna(row["normal_bam"]) and row["normal_bam"] != "":
            bam_path = row["normal_bam"]
            # Only add if we haven't seen this BAM path before
            if bam_path not in seen_bams:
                # Create a unique identifier based on BAM filename
                bam_basename = os.path.basename(bam_path).replace('.bam', '')
                normal_id = f"{bam_basename}_normal"
                # Handle potential filename conflicts by adding counter
                counter = 1
                original_normal_id = normal_id
                while normal_id in normal_bams:
                    normal_id = f"{original_normal_id}_{counter}"
                    counter += 1
                normal_bams[normal_id] = bam_path
                seen_bams.add(bam_path)
    
    return normal_bams

def get_normal_sample_ids():
    """Get list of normal sample identifiers for PoN"""
    return list(get_normal_bams_for_pon().keys())

def get_purecn_candidates():
    """Get samples that can use PureCN (only if PureCN is not skipped)"""
    if config.get("skip_purecn", False):
        return []
    return SAMPLES[(SAMPLES.analysis_type == "TvsN") & (SAMPLES.vcf.notna()) & (SAMPLES.vcf != "")].index.tolist()

def get_vcf_for_call(wildcards):
    try:
        sample_vcf = SAMPLES.loc[wildcards.sample_id, "vcf"]
        # Ensure we have a scalar value, not a Series
        if isinstance(sample_vcf, pd.Series):
            if sample_vcf.empty:
                return []
            sample_vcf = sample_vcf.iloc[0]
        return [sample_vcf] if pd.notna(sample_vcf) and sample_vcf != "" else []
    except KeyError:
        # Sample ID not found in the dataframe
        return []

def get_cnvkit_basename(sample_id):
    """Get the basename that CNVkit will use for output files based on BAM filename"""
    bam_path = SAMPLES.loc[sample_id, "tumor_bam"]
    return os.path.splitext(os.path.basename(bam_path))[0]

# --- Main Workflow Target ---
rule all:
    input:
        expand(f"{config['dirs']['final_vcfs']}/{{sample_id}}.cnv.vcf.gz", sample_id=get_tumor_samples()),
        expand(f"{config['dirs']['plots']}/scatter_genome/{{sample_id}}.genome.pdf", sample_id=get_tumor_samples()),
        expand(f"{config['dirs']['plots']}/diagram/{{sample_id}}.diagram.pdf", sample_id=get_tumor_samples()),
        expand(f"{config['dirs']['plots']}/scatter_gene/{{sample_id}}/{{gene}}.pdf", sample_id=get_tumor_samples(), gene=config["genes_of_interest"]),
        f"{config['dirs']['plots']}/heatmap_cohort.pdf",
        "results/pipeline_summary.csv"

# --- Include Modular Rule Files ---
include: "rules/01_purity.smk"
include: "rules/02_pon.smk"
include: "rules/03_calling.smk"
include: "rules/04_plotting.smk"

# --- Final Summary Report ---
rule create_summary_report:
    input:
        purity_files=expand(f"{config['dirs']['purity_values']}/{{sample_id}}.purity.txt", sample_id=get_tumor_samples()),
        normal_metrics=f"{config['dirs']['pon_creation']}/qc/all_normals.metrics.txt"
    output:
        "results/pipeline_summary.csv"
    log:
        f"{config['dirs']['logs']}/create_summary_report/log.txt"
    run:
        # Purity Summary
        purity_data = []
        for f in input.purity_files:
            sample_id = os.path.basename(f).replace(".purity.txt", "")
            with open(f, 'r') as handle:
                purity = handle.read().strip()
            purity_data.append({"sample_id": sample_id, "final_purity": purity})
        purity_df = pd.DataFrame(purity_data)

        # Normal QC Summary
        try:
            normals_df = pd.read_csv(input.normal_metrics, sep="\t")
            normals_df = normals_df[['sample', 'bivar']].rename(columns={'sample': 'normal_id'})
            normals_df['normal_id'] = normals_df['normal_id'].apply(os.path.basename).str.replace(".targetcoverage.cnn", "")
        except pd.errors.EmptyDataError:
            normals_df = pd.DataFrame(columns=['normal_id', 'bivar']) # Handle empty metrics case
        
        # Merge and save
        if not purity_df.empty:
            purity_df.to_csv(output[0], index=False)
        if not normals_df.empty:
            normals_df.to_csv("results/normal_qc_summary.csv", index=False)
        
        # A more combined report could be built here if desired
        print("Generated summary reports.")
