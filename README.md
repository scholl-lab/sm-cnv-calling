# CNVkit and PureCN Analysis Pipeline

A production-ready, modular Snakemake workflow for calling Copy Number Variations (CNVs) from tumor exome sequencing data. This pipeline integrates **CNVkit** for core coverage analysis and segmentation with **PureCN** for sophisticated tumor purity and ploidy estimation, designed specifically for FFPE tumor samples.

## Features

- **🏗️ Modular Architecture**: Organized into focused modules for maintainability and flexibility
- **🔄 End-to-End Automation**: From BAM files to final VCF calls with comprehensive QC
- **🧬 FFPE-Optimized**: Implements best practices for FFPE data with pooled Panel of Normals and contamination filtering
- **🎯 Hybrid Purity Estimation**: PureCN for Tumor/Normal pairs, fallback estimates for Tumor-only samples
- **📊 Rich Visualizations**: Automated generation of scatter plots, diagrams, heatmaps, and gene-focused plots
- **🐍 Reproducible**: Self-contained Conda environments ensure portability across systems
- **⚡ HPC-Ready**: Optimized for cluster execution with Slurm integration and resource management

## Architecture Overview

The pipeline follows a modular design with five core components:

```
cnv_pipeline.smk (Main Controller)
├── rules/01_purity.smk     # Purity estimation with PureCN
├── rules/02_pon.smk        # Panel of Normals generation
├── rules/03_calling.smk    # CNV calling and export
└── rules/04_plotting.smk   # Comprehensive visualization
```

### Pipeline Workflow

1. **🎯 Target Preparation** (`02_pon.smk`)
   - Annotate targets BED file with gene information from refFlat
   - Split large regions for optimal bin sizes
   - Create antitarget regions for off-target coverage

2. **🧪 Purity Estimation** (`01_purity.smk`)
   - PureCN analysis for `TvsN` samples with VCF files
   - Fallback to user-supplied estimates for other cases
   - Consolidated purity values across all samples

3. **📋 Panel of Normals Generation** (`02_pon.smk`)
   - Parallel CNVkit target and antitarget coverage calculation for all normal samples
   - Automated QC filtering based on noise metrics
   - Pooled reference creation from clean normals

4. **🎯 CNV Calling** (`03_calling.smk`)
   - Tumor sample processing with pooled reference
   - Contamination region exclusion
   - Absolute copy number calling with purity integration
   - VCF export for downstream analysis

5. **📈 Visualization** (`04_plotting.smk`)
   - Scatter plots showing copy number across chromosomes
   - Diagram plots for copy number segments
   - Heatmaps for cross-sample comparison
   - Gene-focused plots for regions of interest

## Quick Start

### 1. Installation

```bash
git clone [repository-url]
cd sm-cnv-calling
```

**Prerequisites**: [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda installed

#### A. Create Conda Environments

The pipeline uses two main conda environments that need to be created before running. You can either let Snakemake create them automatically or set them up manually:

**Option 1: Automatic Setup (Recommended)**
```bash
# Snakemake will automatically create the environments when first run
snakemake -s scripts/snakemake/cnv_pipeline.smk --use-conda --conda-prefix conda_envs --cores 1 --dry-run
```

**Option 2: Manual Setup**
```bash
# Create the conda environments manually with specific names
conda env create -f scripts/snakemake/envs/cnvkit.yaml -n cnvkit
conda env create -f scripts/snakemake/envs/purecn.yaml -n purecn

# Alternatively, create them in the conda_envs directory as Snakemake expects
mkdir -p conda_envs
conda env create -f scripts/snakemake/envs/cnvkit.yaml -p conda_envs/cnvkit
conda env create -f scripts/snakemake/envs/purecn.yaml -p conda_envs/purecn
```

The pipeline will automatically use these environments for different steps:
- **`cnvkit`**: Used for CNVkit operations, coverage analysis, and plotting
- **`purecn`**: Used for PureCN-based purity estimation

### 2. Configuration

#### B. Download Standard Reference Files (Optional)

This repository includes a helper script to download standard accessory files (`access*.bed`, `refFlat.txt`) for common genome builds.

```bash
# Create a directory for your references
mkdir -p references

# Activate the cnvkit conda environment (if created manually)
# conda activate cnvkit
# Or if using conda_envs directory:
# conda activate conda_envs/cnvkit

# Download files for hg19
python scripts/snakemake/helpers/download_references.py --genome hg19 --output-dir references/
```

This will download `access-5k-mappable.hg19.bed` and `refFlat.txt` into the `references/` directory.

#### C. Configure the Pipeline

Edit `scripts/snakemake/config.yaml` to specify your reference files. If you used the downloader script, these paths would be `references/access-5k-mappable.hg19.bed`, `references/refFlat.txt`, etc.:

```yaml
reference_genome: "/path/to/your/ucsc.hg19.fasta"
targets_bed: "/path/to/your/targets.bed"
access_bed: "references/access-5k-mappable.hg19.bed"
blacklist: "/path/to/your/genomic_blacklist.bed"
annotate_refFlat: "references/refFlat.txt"  # Optional: comment out to disable gene annotation
purecn_normal_panel_vcf: "/path/to/your/normal_panel.vcf.gz"

# Optional: Skip PureCN analysis entirely (default: false)
skip_purecn: false
```

**Skipping PureCN**: If you encounter issues with PureCN dependencies or want to run only CNVkit analysis, set `skip_purecn: true` in the config. This will:
- Skip all PureCN-related rules and dependencies
- Use only fallback purity values from the sample sheet
- Still run the complete CNVkit workflow for CNV calling

#### D. Prepare Sample Sheet

Create your `scripts/snakemake/samples.tsv` with sample information (see [Input Formats](#input-formats) below).

### 3. Execution

The conda environments will be automatically created by Snakemake when using the `--use-conda` flag. Execute the pipeline:

```bash
# Run the complete pipeline (environments created automatically)
snakemake -s scripts/snakemake/cnv_pipeline.smk --use-conda --conda-prefix conda_envs --cores 8

# For cluster execution (SLURM)
sbatch scripts/run_cnv_pipeline.sh
```

### 4. Using the Reference Downloader

For additional reference files (like hg38 or updated files), you can use the included downloader:

```bash
# Download hg38 files
python scripts/snakemake/helpers/download_references.py --genome hg38 --output-dir references/

# Download with custom config path
python scripts/snakemake/helpers/download_references.py --genome hg19 --output-dir references/ --config my_custom_download_config.yaml
```

## Configuration Reference

### Core Configuration (`config.yaml`)

#### Directory Structure
All output paths are centralized under the `dirs` section:

```yaml
dirs:
  purecn_setup: "results/01_purecn_setup"      # PureCN interval files
  purecn_runs: "results/02_purecn_runs"        # PureCN analysis outputs
  purity_values: "results/03_purity_values"    # Consolidated purity estimates
  pon_creation: "results/04_pon_creation"      # Panel of Normals files
  cnvkit_runs: "results/05_cnvkit_runs"        # CNVkit coverage/segmentation
  final_calls: "results/06_final_calls"        # Final copy number calls
  final_vcfs: "results/07_final_vcfs"          # Exported VCF files
  plots: "results/08_plots"                    # All visualization outputs
  logs: "logs"                                 # Pipeline execution logs
```

#### Reference Files (Required)
Update these paths to match your system:

| Parameter | Description | Format |
|-----------|-------------|---------|
| `reference_genome` | Reference FASTA file | `.fasta`, `.fa` |
| `targets_bed` | Target regions for exome capture | `.bed` |
| `access_bed` | Mappable regions (e.g., UCSC access-5k) | `.bed` |
| `blacklist` | Contaminated regions to exclude | `.bed` |
| `annotate_refFlat` | Gene annotations for CNVkit (optional) | `.txt` |
| `purecn_normal_panel_vcf` | PureCN normal panel VCF | `.vcf.gz` |

#### Conda Environment Configuration

The pipeline uses two conda environments defined in the `conda_envs` section:

```yaml
conda_envs:
  cnvkit: "scripts/snakemake/envs/cnvkit.yaml"    # CNVkit operations and plotting
  purecn: "scripts/snakemake/envs/purecn.yaml"    # PureCN purity estimation
```

When using `--conda-prefix conda_envs`, Snakemake will create environments in:
- `conda_envs/cnvkit/` - Contains cnvkit, pandas, requests, tqdm, and related tools
- `conda_envs/purecn/` - Contains R, PureCN, and related Bioconductor packages

#### Tool Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `purecn_genome` | `"hg19"` | Reference genome build |
| `purecn_assay_name` | `"MyTumorAssay"` | Assay identifier for PureCN |
| `pon_qc_metric_threshold` | `0.3` | Maximum noise threshold for normals |
| `segment_threshold` | `1.0e-4` | P-value threshold for segmentation |

#### Visualization Parameters

```yaml
genes_of_interest:
  - "TP53"
  - "MYC"
  - "EGFR"
  - "ERBB2"
```

Specify genes for focused plotting. Gene names should match your reference annotation.

#### Resource Allocation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `default_threads` | `4` | Default CPU cores per job |
| `default_mem_mb` | `16000` | Default memory (MB) per job |
| `cnvkit_batch_mem_mb` | `24000` | Memory for CNVkit batch processing |
| `cnvkit_call_mem_mb` | `8000` | Memory for CNVkit call processing with VCFs |
| `cnvkit_export_mem_mb` | `4000` | Memory for CNVkit VCF export |
| `cnvkit_plot_mem_mb` | `6000` | Memory for CNVkit individual plot operations |
| `cnvkit_heatmap_mem_mb` | `24000` | Memory for CNVkit cohort heatmap generation (increase to 32000+ for >50 samples) |

## Input Formats

### Sample Sheet (`samples.tsv`)

The pipeline requires a tab-separated file with the following columns:

| Column | Required | Description | Example |
|--------|----------|-------------|---------|
| `sample_id` | ✅ | Unique sample identifier | `PatientA` |
| `analysis_type` | ✅ | Analysis mode (see below) | `TvsN` |
| `tumor_bam` | Conditional | Path to tumor BAM file | `/path/to/tumor.bam` |
| `normal_bam` | Conditional | Path to normal BAM file | `/path/to/normal.bam` |
| `vcf` | Optional | Path to variant VCF file | `/path/to/variants.vcf.gz` |
| `tumor_sample_id_vcf` | Optional | Tumor sample ID in VCF for B-allele frequency calculation | `PatientA-T` |
| `normal_sample_id_vcf` | Optional | Normal sample ID in VCF for B-allele frequency calculation | `PatientA-N` |
| `fallback_purity` | ✅ | Backup purity estimate (0-1) | `0.35` |

#### Analysis Types

| Type | Description | Required Files | Use Case |
|------|-------------|----------------|----------|
| `TvsN` | Tumor vs Normal | `tumor_bam`, `normal_bam`, `vcf` (optional) | Paired analysis with germline control; normal BAM used for Panel of Normals |
| `To` | Tumor only | `tumor_bam`, `vcf` (optional) | Single sample analysis without matched normal |

#### Example Sample Sheet

```tsv
sample_id	analysis_type	tumor_bam	normal_bam	vcf	tumor_sample_id_vcf	normal_sample_id_vcf	fallback_purity
PatientA_T1	TvsN	/data/bams/PatientA_T1.bam	/data/bams/PatientA_N.bam	/data/vcfs/PatientA_T1.vcf.gz	PatientA-T	PatientA-N	0.35
PatientA_T2	TvsN	/data/bams/PatientA_T2.bam	/data/bams/PatientA_N.bam	/data/vcfs/PatientA_T2.vcf.gz	PatientA-T	PatientA-N	0.40
PatientB	TvsN	/data/bams/PatientB_T.bam	/data/bams/PatientB_N.bam	/data/vcfs/PatientB.vcf.gz	PatientB-T	PatientB-N	0.30
PatientC	To	/data/bams/PatientC_T.bam		/data/vcfs/PatientC.vcf.gz	PatientC-T		0.45
```

**Note**: 
- The normal BAM files from `TvsN` samples are automatically used to construct the Panel of Normals (PoN) for improved CNV calling accuracy.
- **Multiple tumors per patient**: If the same normal BAM is used for multiple tumor samples (e.g., `PatientA_T1` and `PatientA_T2` both use `PatientA_N.bam`), the pipeline automatically deduplicates normal BAMs to avoid processing the same file multiple times in the PoN.
- **VCF Sample IDs**: For samples with VCF files, add `tumor_sample_id_vcf` and `normal_sample_id_vcf` columns to specify the exact sample identifiers in the VCF file. This enables proper B-allele frequency calculation from germline variants and accurate VCF parsing when sample names in the VCF differ from the sample sheet IDs.

### BAM File Requirements

- **Coordinate-sorted** and **indexed** (`.bai` files required)
- **Aligned to the same reference** specified in `config.yaml`
- **Duplicate-marked** (recommended but not required)
- **Base quality scores recalibrated** (recommended for optimal results)

### VCF File Requirements

- **Compressed and indexed** (`.gz` + `.tbi` files)
- **Same reference coordinate system** as BAM files
- **Contains germline variants** for optimal PureCN analysis
- **Quality-filtered variants** recommended (PASS filter)

## Output Structure

The pipeline generates organized outputs across multiple directories:

### 📊 Visualization Outputs (`results/08_plots/`)

- **`scatter_genome/{sample}.genome.pdf`**: Copy number scatter plots across chromosomes
- **`diagram/{sample}.diagram.pdf`**: Segmented copy number diagrams  
- **`scatter_gene/{sample}/{gene}.pdf`**: Gene-focused plots for regions of interest, organized by sample
- **`heatmap_cohort.pdf`**: Cross-sample copy number heatmap

**Note**: Gene-focused plots are organized into sample-specific subdirectories for better organization. For example, with samples `PatientA_T1`, `PatientB` and genes `TP53`, `MYC`, the structure would be:
```
results/08_plots/scatter_gene/
├── PatientA_T1/
│   ├── TP53.pdf
│   └── MYC.pdf
└── PatientB/
    ├── TP53.pdf
    └── MYC.pdf
```

### 📄 Analysis Results

- **`results/07_final_vcfs/`**: Final VCF files with copy number calls
- **`results/06_final_calls/`**: CNVkit native format calls (`.call.cns`)
- **`results/03_purity_values/`**: Consolidated purity estimates
- **`results/04_pon_creation/pooled_reference.cnn`**: Panel of Normals reference

### 🔍 Intermediate Files

- **`results/01_purecn_setup/`**: PureCN interval and setup files
- **`results/02_purecn_runs/`**: Individual PureCN analysis results
- **`results/05_cnvkit_runs/`**: CNVkit coverage and segmentation files

## Advanced Usage

### Running Specific Modules

Execute only purity estimation:
```bash
snakemake -s scripts/snakemake/cnv_pipeline.smk --use-conda \
  purity_consolidation \
  --cores 8
```

Generate only visualizations:
```bash
snakemake -s scripts/snakemake/cnv_pipeline.smk --use-conda \
  all_plots \
  --cores 4
```

### Custom Resource Allocation

Override default resources in your execution:
```bash
snakemake -s scripts/snakemake/cnv_pipeline.smk --use-conda \
  --set-resources cnvkit_batch:mem_mb=32000 \
  --cores 16
```

### Cluster Execution

For SLURM clusters, use the provided execution script:
```bash
# Edit cluster parameters in run_cnv_pipeline.sh
./scripts/run_cnv_pipeline.sh
```

Or create a custom SLURM profile following [Snakemake cluster documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

### Debugging and Troubleshooting

#### Check Pipeline Status
```bash
snakemake -s scripts/snakemake/cnv_pipeline.smk --summary
```

#### Validate Configuration
```bash
snakemake -s scripts/snakemake/cnv_pipeline.smk --lint
```

#### Generate Pipeline Visualization
```bash
snakemake -s scripts/snakemake/cnv_pipeline.smk --dag | dot -Tpng > pipeline_dag.png
```

## Quality Control

The pipeline includes several QC checkpoints:

1. **Reference File Validation**: Automatic verification of all required reference files
2. **Normal Sample QC**: Automated filtering of noisy normal samples using CNVkit metrics
3. **Purity Validation**: Consolidated purity estimates with fallback handling
4. **Output Validation**: File existence and format checks throughout the pipeline

## Performance Considerations

### Memory Requirements

- **Minimum**: 16 GB RAM for small exomes
- **Recommended**: 32 GB RAM for large panels or whole exomes
- **PureCN**: May require up to 24 GB for complex samples

### Runtime Estimates

| Step | Single Sample | 10 Samples | 50 Samples |
|------|---------------|------------|------------|
| Purity Estimation | 30-60 min | 2-4 hours | 6-12 hours |
| Panel of Normals | N/A | 1-2 hours | 3-6 hours |
| CNV Calling | 15-30 min | 1-3 hours | 4-8 hours |
| Visualization | 5-10 min | 30-60 min | 2-4 hours |

*Estimates based on typical exome data (~50 MB BAM files) running on modern HPC infrastructure*

## Citation

If you use this pipeline in your research, please cite:

- **CNVkit**: Talevich et al. (2016). CNVkit: Genome-wide copy number detection and visualization from targeted sequencing. *PLOS Computational Biology*
- **PureCN**: Riester et al. (2016). PureCN: copy number calling and SNV classification using targeted short read sequencing. *Source Code for Biology and Medicine*

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all existing tests pass
5. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For questions and support:

- **Issues**: Use GitHub Issues for bug reports and feature requests
- **Discussions**: Use GitHub Discussions for usage questions and community support
- **Documentation**: Refer to individual tool documentation for [CNVkit](https://cnvkit.readthedocs.io/) and [PureCN](https://bioconductor.org/packages/release/bioc/html/PureCN.html)

Submit the main pipeline job using the provided shell script:
```bash
# Submit with default of 50 parallel jobs
sbatch scripts/run_cnv_pipeline.sh

# Submit with a custom number of parallel jobs
sbatch scripts/run_cnv_pipeline.sh 100
```

The script will handle temporary directories and logging. All final results will be placed in the `results/` directory, following the structure outlined in the `Snakefile`.

