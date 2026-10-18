# Deprecated Legacy Workflow Files

This directory archives the legacy Snakemake 7 `scripts/snakemake/` pipeline and `scripts/run_cnv_pipeline.sh` launcher.

## Migration to Snakemake 8+ Layout

- **Entry point**: `workflow/Snakefile` (replaces `scripts/snakemake/cnv_pipeline.smk`)
- **Rules & helpers**: `workflow/rules/{common.smk,helpers.py,cnvkit.smk,purecn.smk,qc.smk}` (replaces `scripts/snakemake/rules/*.smk` and `scripts/snakemake/helpers/*.py`)
- **Conda environments**: `workflow/envs/{cnvkit.yaml,purecn.yaml}` (replaces `scripts/snakemake/envs/*.yaml`)
- **Configuration & schemas**: `config/config.yaml`, `config/samples.tsv`, and `workflow/schemas/{config.schema.yaml,samples.schema.yaml}`
- **SLURM launcher & profiles**: `scripts/run_snakemake.sh` and `profiles/{default,bih,charite,local}/config.yaml` (replaces `scripts/run_cnv_pipeline.sh`)
