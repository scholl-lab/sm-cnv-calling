# CNVkit and PureCN Analysis Pipeline

This repository contains a robust, reproducible Snakemake workflow for calling Copy Number Variations (CNVs) from tumor exome sequencing data. It integrates **CNVkit** for core coverage analysis and segmentation with **PureCN** for sophisticated tumor purity and ploidy estimation.

## Features

- **End-to-End Automation**: From BAM files to final VCF calls.
- **Best Practices**: Implements best practices for FFPE data, including a pooled Panel of Normals (PoN) with automated QC, and blacklisting of contaminated regions.
- **Hybrid Purity Estimation**: Uses PureCN for Tumor/Normal pairs with VCFs and provides a fallback to user-supplied estimates for Tumor-only cases.
- **Reproducibility**: Uses Conda environments to manage all software dependencies, ensuring the pipeline is portable and results are reproducible.
- **Scalability**: Designed for HPC clusters using a Slurm submission script and a Snakemake profile.

## Pipeline Overview

1.  **Purity Estimation**:
    - For `TvsN` samples with a VCF, PureCN is run to estimate purity.
    - For all other tumor samples (`To` or `TvsN` without a VCF), a fallback purity from the `samples.tsv` is used.
2.  **Panel of Normals (PoN) Generation**:
    - All samples marked as `Normal` are processed with CNVkit.
    - `cnvkit.py metrics` is used to assess the noise profile of each normal.
    - A helper script filters out high-noise normals based on a configurable threshold.
    - A final pooled reference (`pooled_reference.cnn`) is built from the clean set of normals.
3.  **CNV Calling**:
    - Each tumor sample is processed with `cnvkit.py batch` using the robust pooled reference.
    - Contaminated regions specified in the config are excluded.
    - `cnvkit.py call` is used to determine absolute integer copy numbers, using the final purity estimate and B-allele frequencies from the VCF.
4.  **Final Export**:
    - The final calls are exported to a compressed VCF format for downstream analysis.

## Getting Started

### 1. Installation

Clone this repository:
```bash
git clone [repository-url]
cd cnv_pipeline
```

Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda installed. The pipeline will automatically create the necessary environments.

### 2. Configuration

Modify the following files in `scripts/snakemake/`:

1.  **`config.yaml`**: Update all paths to reference files (`reference_genome`, `targets_bed`, etc.) and adjust any tool parameters as needed.
2.  **`samples.tsv`**: This is the most important file. Populate it with your sample information. Provide a unique `sample_id` for each tumor, the `analysis_type`, paths to `tumor_bam`, `normal_bam`, and `vcf` files, and a `fallback_purity`.

### 3. Execution

The pipeline is designed to be run on an HPC cluster with a Slurm workload manager.

Submit the main pipeline job using the provided shell script:
```bash
# Submit with default of 50 parallel jobs
sbatch scripts/run_cnv_pipeline.sh

# Submit with a custom number of parallel jobs
sbatch scripts/run_cnv_pipeline.sh 100
```

The script will handle temporary directories and logging. All final results will be placed in the `results/` directory, following the structure outlined in the `Snakefile`.

