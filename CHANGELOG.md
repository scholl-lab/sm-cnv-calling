# Changelog

All notable changes to `sm-cnv-calling` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **PoN noise QC used raw coverage, so no normal ever passed (flat reference):** `cnvkit_metrics_normals` ran `cnvkit.py metrics` on raw `*.targetcoverage.cnn` files, whose spread is large for every normal (apa-sequencing run v3: bivar 0.46-0.55 for all 20 normals vs. `pon_bivar_threshold` 0.3). `identify_clean_normals` selected 0 normals and `cnvkit_reference_pooled` silently built a flat reference, giving noisy segmentation, impossible chrX copy numbers and poor PureCN fits. **All earlier runs that had normals used a flat reference; re-run `cnvkit_reference_pooled` and everything downstream.** New rules `cnvkit_reference_prelim` (reference from all normals, `pon_qc/prelim_reference.cnn`) and `cnvkit_fix_normal` (`pon_qc/{normal}.cnr`) normalise each normal first (self-included, so it ranks outliers rather than measuring absolute noise); `cnvkit_metrics_normals` now runs on those `.cnr` files and `identify_clean_normals` also strips `.cnr`.
### Changed
- `cnvkit_reference_pooled` builds a flat reference only when the cohort has no normal BAMs. If normals exist but none passes `pon_bivar_threshold`, it now fails with an error instead of silently falling back; the log states how many normals passed.

## [0.1.3] - 2026-10-23

### Fixed
- **`cnvkit_segment` corrupt segments with VCF BAF (CNVkit 0.9.14):** `segment -m cbs -v VCF` concatenates the per-segment BAF results without `ignore_index`, so duplicate row labels make `transfer_fields` stretch many sub-segments to whole-arm coordinates and broadcast BAF across rows (APA15: 82 overlapping segments, duplicate coordinates; `cnvkit.py heatmap` crashed with "cannot handle a non-unique multi-index"). The rule now uses a new env `workflow/envs/cnvkit_segment.yaml` (same as `cnvkit.yaml`) whose `cnvkit_segment.post-deploy.sh` installs etal/cnvkit commit `7539ce3` (0.9.15.dev0, includes the upstream fix etal/cnvkit#1125) with `pip install --no-deps` and fails if the version differs. Verified on APA15: 0 overlaps, 0 duplicates, same probe total, heatmap OK. PureCN uses the CNVkit `.seg` export, so purity/ploidy fits made on 0.9.14 segments were affected; re-run `cnvkit_segment` and downstream. Drop the env once a CNVkit release newer than 0.9.14 contains #1125.
- **`purecn_run` load balancing:** with `--cores N`, `PureCN.R` builds `MulticoreParam(workers = N)` with BiocParallel's default `tasks = 0`, which pre-splits the ~17 purity/ploidy local optima into N fixed chunks; one worker ran alone for over 19 h while the others idled (seen in apa-sequencing run v3, APA55). When `threads > 1` the rule now runs a copy of `PureCN.R` patched to `tasks = 1000L` (capped at the number of optima), so each optimum is its own task. The rule fails loudly if the expected line is not found in the installed PureCN.

## [0.1.2] - 2026-10-21

### Fixed
- **PureCN on multi-region units:** each multi-region unit (`{sample}__{bam}`) got the joint multi-sample Mutect2 VCF (e.g. two tumour regions plus the normal), and `PureCN.R` stopped with `VCF contains 3 samples. Should be either one or two`. The fault-tolerant fallback then wrote empty `.csv`/`.rds` with status `failed`, so the failure was silent. The new rule `purecn_subset_vcf` (bcftools, CNVkit env) subsets the VCF to the unit's `tumor_vcf_id` (+ `normal_vcf_id`) and PureCN uses the subset. Single-region units keep their VCF input unchanged, so existing PureCN results are not invalidated. Units carry a `multi_region` flag.
### Changed
- **`purecn_run` uses its threads:** it passes `--cores {threads}` to `PureCN.R` (BiocParallel MulticoreParam), so the purity/ploidy local optima are fitted in parallel. Set the thread count in the profile (`set-threads: purecn_run`); the default stays 1. Serial fitting took about 10 min per candidate on FFPE exomes.
- **`unit` wildcard constraint (#15):** the constraint is now `(?:(?!\.segmetrics)[^/])+`, which checks only the wildcard's own characters, so `cnvkit_segment` no longer tries `unit=<id>.segmetrics`. This removes one `InputFunctionException ... KeyError` per unit from every DAG build.

## [0.1.0] - 2026-10-19

### Added
- **Snakemake 8+ Workflow Layout & Schemas (#4)**:
  - Migrated workflow to `workflow/Snakefile`, `workflow/rules/{common,cnvkit,purecn,qc}.smk`, and `workflow/rules/helpers.py` with `min_version("8.0")`.
  - Added JSON Schemas `workflow/schemas/config.schema.yaml` and `workflow/schemas/samples.schema.yaml` with `GRCh38`/`hg38` defaults.
  - Added pinned Bioconda environments `workflow/envs/cnvkit.yaml` (`cnvkit=0.9.14`, `bedtools=2.31.1`, `samtools=1.21`, `bcftools=1.21`, `htslib=1.21`) and `workflow/envs/purecn.yaml` (`bioconductor-purecn=2.16.1`, `r-pscbs=0.68.0`).
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
