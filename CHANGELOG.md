# Changelog

All notable changes to `sm-cnv-calling` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-10-19

### Added
- **Snakemake 8+ Workflow Layout & Schemas (#4)**:
  - Migrated workflow to `workflow/Snakefile`, `workflow/rules/{common,cnvkit,purecn,qc}.smk`, and `workflow/rules/helpers.py` with `min_version("8.0")`.
  - Added JSON Schemas `workflow/schemas/config.schema.yaml` and `workflow/schemas/samples.schema.yaml` with `GRCh38`/`hg38` defaults.
  - Added pinned Bioconda environments `workflow/envs/cnvkit.yaml` (`cnvkit=0.9.12`, `bedtools=2.31.1`, `samtools=1.21`, `bcftools=1.21`, `htslib=1.21`) and `workflow/envs/purecn.yaml` (`bioconductor-purecn=2.12.0`, `bioconductor-dnacopy=1.78.0`, `r-pscbs=0.67.0`).
  - Added `scripts/generate_config.py` for automated sample sheet and config generation from `sm-alignment` BAM directories.
  - Archived legacy `scripts/snakemake` and `scripts/run_cnv_pipeline.sh` under `deprecated/`.
- **Shared Sample Sheet & Multi-Region / 97-BAM Cohort Support (#5)**:
  - Aligned `config/samples.tsv` schema and BAM/VCF path resolution with `sm-calling` (`tumor_normal`, `tumor_only`, `germline`, plus legacy `TvsN`/`To`).
  - Added expansion of `;`-delimited multi-region `tumor_bam` entries into `{sample}__{bam_basename}` units with per-unit VCF sample ID resolution (`-i`/`-n`).
- **Accessible Genome & Blacklist Subtraction (#1)**:
  - Added automatic `cnvkit.py access` generation when `ref.access_bed` is omitted and `bedtools subtract` blacklist removal for both target and access BEDs.
- **Target Bin Splitting, Gene Annotation & Autobin (#2)**:
  - Updated `cnvkit_target` to always run `cnvkit.py target --split --short-names` with optional `--annotate {ref_flat}` and `--avg-size`, plus optional `cnvkit.autobin` support.
- **Modular Per-BAM Coverage, Allele-Aware CBS & Bintest (#3)**:
  - Replaced monolithic `cnvkit_batch_tumor` and `cnvkit_for_purecn` with per-unique-BAM `cnvkit_coverage_target` and `cnvkit_coverage_antitarget`, per-unit `cnvkit_fix`, BAF-aware `cnvkit_segment` (`-m cbs -t 1e-4 --drop-low-coverage`), and `cnvkit_bintest`.
- **Chromosomal Sex Inference, Segmetrics CIs & Standard Exports (#8)**:
  - Added `cnvkit_sex` inference and sex conflict checking (`-y` and `-x` propagated across reference, calling, diagrams, and exports).
  - Added `cnvkit_segmetrics` (`--ci --pi --sem --iqr`), `cnvkit_call --filter ci`, `cnvkit_metrics_tumor`, `cnvkit_export_seg`, `cnvkit_export_bed`, `cnvkit_export_vcf` (bgzipped + tabix-indexed), and guarded `cnvkit_scatter_gene`.
- **Snakemake 8+ SLURM Profiles & Launcher (#9)**:
  - Added `profiles/{default,bih,charite,local}/config.yaml` and cluster-detecting `scripts/run_snakemake.sh`.
- **Automated Test Suite, Linting & CI (#10)**:
  - Added `pyproject.toml` (`version = "0.1.0"`), `Makefile`, `.github/workflows/ci.yml`, and `tests/{conftest,test_schema_validation,test_helpers,test_dryrun}.py`.

### Fixed
- **Panel of Normals (PoN) QC & Flat Reference Fallback (#6)**:
  - Fixed normal sample-to-`.cnn` path mapping in `identify_clean_normals` to output verified paths for both `.targetcoverage.cnn` and `.antitargetcoverage.cnn`.
  - Fixed empty-input crash in `cnvkit_metrics_normals` and added automatic flat reference fallback when 0 clean normals are available.
- **PureCN Integration & Purity/Ploidy Propagation (#7)**:
  - Fixed PureCN integration to support both `cnvkit` mode (PoN-normalized `.cnr` + `.seg`) and `native` mode (`IntervalFile.R --offtarget`, `Coverage.R`, `NormalDB.R --coveragefiles`) for `tumor_normal` and `tumor_only` samples.
  - Propagated both `purity` and `ploidy` from PureCN to `cnvkit_call` and replaced fabricated `fallback_purity` with safe fallback to `cnvkit.py call -m threshold`.
