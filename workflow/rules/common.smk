"""Shared configuration shortcuts, metadata accessors, and helper wrappers."""

import os
from pathlib import Path

from rules.helpers import (
    build_cnvkit_vcf_args as _build_cnvkit_vcf_args_impl,
    get_bam_path_map as _get_bam_path_map_impl,
    get_cnv_unit_ids as _get_cnv_unit_ids_impl,
    get_cnv_units as _get_cnv_units_impl,
    get_germline_samples as _get_germline_samples_impl,
    get_normal_bam_basenames as _get_normal_bam_basenames_impl,
    get_samples as _get_samples_impl,
    get_tumor_bam_basenames as _get_tumor_bam_basenames_impl,
    get_tumor_samples as _get_tumor_samples_impl,
    get_unique_bam_basenames as _get_unique_bam_basenames_impl,
    has_matched_normal as _has_matched_normal_impl,
    resolve_sample_vcf as _resolve_sample_vcf_impl,
    resolve_sample_vcf_stats as _resolve_sample_vcf_stats_impl,
)

try:
    RULES_DIR = str(_WORKFLOW_DIR / "rules")
except NameError:
    RULES_DIR = str(Path(__file__).resolve().parent)


# -- Reference configuration -------------------------------------------------
_ref_cfg = config["ref"]
REF = _ref_cfg["genome"]
REF_BUILD = _ref_cfg["build"]
TARGETS_BED = _ref_cfg["targets_bed"]
BAITS_BED = _ref_cfg.get("baits_bed", "") or TARGETS_BED
REF_ACCESS_BED = _ref_cfg.get("access_bed", "")
REF_BLACKLIST_BED = _ref_cfg.get("blacklist_bed", "")
REF_FLAT = _ref_cfg.get("ref_flat", "")
REF_MAPPABILITY_BW = _ref_cfg.get("mappability_bw", "")
REF_REPTIMING_BW = _ref_cfg.get("reptiming_bw", "")


# -- Paths & BAM settings ----------------------------------------------------
_paths_cfg = config["paths"]
BAM_FOLDER = _paths_cfg["bam_folder"]
VCF_FOLDER = _paths_cfg.get("vcf_folder", "")
OUTPUT_DIR = _paths_cfg["output_folder"]
LOG_SUBDIR = _paths_cfg.get("log_subdir", "logs")
LOG_DIR = os.path.join(OUTPUT_DIR, LOG_SUBDIR)

BAM_EXT = config["bam"]["file_extension"]


# -- CNVkit configuration ----------------------------------------------------
_cnvkit_cfg = config.get("cnvkit", {})
_cnvkit_params = config.get("params", {}).get("cnvkit", {})

CNVKIT_REFERENCE_CNN = _cnvkit_cfg.get("reference_cnn", "")
CNVKIT_AUTOBIN = bool(_cnvkit_cfg.get("autobin", False))
CNVKIT_TARGET_AVG_SIZE = int(_cnvkit_cfg.get("target_avg_size", 200))
CNVKIT_ANTITARGET_AVG_SIZE = int(_cnvkit_cfg.get("antitarget_avg_size", 150000))
CNVKIT_ANTITARGET_MIN_SIZE = int(_cnvkit_cfg.get("antitarget_min_size", 500))
CNVKIT_ACCESS_MIN_GAP_SIZE = int(_cnvkit_cfg.get("access_min_gap_size", 5000))
CNVKIT_PON_BIVAR_THRESHOLD = float(_cnvkit_cfg.get("pon_bivar_threshold", 0.3))
CNVKIT_MALE_REFERENCE = bool(_cnvkit_cfg.get("male_reference", False))
CNVKIT_SEGMENT_METHOD = _cnvkit_params.get(
    "segment_method", _cnvkit_cfg.get("segment_method", "cbs")
)
CNVKIT_SEGMENT_THRESHOLD = float(
    _cnvkit_params.get("segment_threshold", _cnvkit_cfg.get("segment_threshold", 1.0e-4))
)
CNVKIT_DROP_LOW_COVERAGE = bool(_cnvkit_cfg.get("drop_low_coverage", True))
CNVKIT_CALL_METHOD = _cnvkit_cfg.get("call_method", "clonal")
CNVKIT_CALL_FILTER = _cnvkit_cfg.get("call_filter", "ci")
CNVKIT_MIN_VARIANT_DEPTH = int(_cnvkit_cfg.get("min_variant_depth", 20))
CNVKIT_GENES = list(_cnvkit_cfg.get("genes_of_interest", ["TP53", "MYC", "EGFR", "ERBB2"]))

CNVKIT_ACCESS_EXTRA = _cnvkit_params.get("access_extra", "")
CNVKIT_TARGET_EXTRA = _cnvkit_params.get("target_extra", "")
CNVKIT_ANTITARGET_EXTRA = _cnvkit_params.get("antitarget_extra", "")
CNVKIT_AUTOBIN_EXTRA = _cnvkit_params.get("autobin_extra", "")
CNVKIT_COVERAGE_EXTRA = _cnvkit_params.get("coverage_extra", "")
CNVKIT_REFERENCE_EXTRA = _cnvkit_params.get("reference_extra", "")
CNVKIT_FIX_EXTRA = _cnvkit_params.get("fix_extra", "")
CNVKIT_SEGMENT_EXTRA = _cnvkit_params.get("segment_extra", "")
CNVKIT_SEGMETRICS_EXTRA = _cnvkit_params.get("segmetrics_extra", "")
CNVKIT_CALL_EXTRA = _cnvkit_params.get("call_extra", "")
CNVKIT_BINTEST_EXTRA = _cnvkit_params.get("bintest_extra", "")
CNVKIT_SCATTER_EXTRA = _cnvkit_params.get("scatter_extra", "")
CNVKIT_DIAGRAM_EXTRA = _cnvkit_params.get("diagram_extra", "")
CNVKIT_HEATMAP_EXTRA = _cnvkit_params.get("heatmap_extra", "")


# -- PureCN configuration ----------------------------------------------------
_purecn_cfg = config.get("purecn", {})
_purecn_params = config.get("params", {}).get("purecn", {})

PURECN_ENABLED = bool(_purecn_cfg.get("enabled", True))
PURECN_MODE = _purecn_cfg.get("mode", "cnvkit")
PURECN_GENOME = _purecn_cfg.get("genome", "hg38")
PURECN_ASSAY = _purecn_cfg.get("assay", "exome")
PURECN_FUN_SEGMENTATION = _purecn_cfg.get("fun_segmentation", "none")
PURECN_NORMAL_PANEL_VCF = _purecn_cfg.get("normal_panel_vcf", "")
PURECN_NORMALDB = _purecn_cfg.get("normaldb", "")
PURECN_MAPPING_BIAS = _purecn_cfg.get("mapping_bias", "")
PURECN_SNP_BLACKLIST = _purecn_cfg.get("snp_blacklist", "")
PURECN_SEED = int(_purecn_cfg.get("seed", 123))
PURECN_POSTOPTIMIZE = bool(_purecn_cfg.get("postoptimize", True))
PURECN_EXTRA = _purecn_params.get("run_extra", _purecn_cfg.get("extra", ""))
PURECN_INTERVALFILE_EXTRA = _purecn_params.get("intervalfile_extra", "")
PURECN_COVERAGE_EXTRA = _purecn_params.get("coverage_extra", "")
PURECN_NORMALDB_EXTRA = _purecn_params.get("normaldb_extra", "")


# -- Output directories ------------------------------------------------------
CNVKIT_DIR = os.path.join(OUTPUT_DIR, "cnvkit")
CNVKIT_REF_DIR = os.path.join(CNVKIT_DIR, "reference")
CNVKIT_COV_DIR = os.path.join(CNVKIT_DIR, "coverage")
CNVKIT_PON_QC_DIR = os.path.join(CNVKIT_DIR, "pon_qc")
CNVKIT_RUNS_DIR = os.path.join(CNVKIT_DIR, "runs")
CNVKIT_CALLS_DIR = os.path.join(CNVKIT_DIR, "calls")
CNVKIT_EXPORTS_DIR = os.path.join(CNVKIT_DIR, "exports")
CNVKIT_PLOTS_DIR = os.path.join(CNVKIT_DIR, "plots")

PURECN_DIR = os.path.join(OUTPUT_DIR, "purecn")
PURECN_REF_DIR = os.path.join(PURECN_DIR, "reference")
PURECN_COV_DIR = os.path.join(PURECN_DIR, "coverage")
PURECN_RUNS_DIR = os.path.join(PURECN_DIR, "runs")
PURITY_DIR = os.path.join(PURECN_DIR, "purity_ploidy")

QC_DIR = os.path.join(OUTPUT_DIR, "qc")


# -- Sample & unit accessors -------------------------------------------------
BAM_PATH_MAP = _get_bam_path_map_impl(samples_df, BAM_FOLDER, BAM_EXT)
CNV_UNITS = _get_cnv_units_impl(samples_df, BAM_EXT)


wildcard_constraints:
    unit=r"(?!.*\.segmetrics$)[^/]+",
    bam_basename=r"[^/]+",


def get_all_samples():
    return _get_samples_impl(samples_df)


def get_tumor_samples():
    return _get_tumor_samples_impl(samples_df)


def get_germline_samples():
    return _get_germline_samples_impl(samples_df)


def get_unique_bams():
    return _get_unique_bam_basenames_impl(samples_df, BAM_EXT)


def get_normal_bams():
    return _get_normal_bam_basenames_impl(samples_df, BAM_EXT)


def get_tumor_bams(sample):
    return _get_tumor_bam_basenames_impl(samples_df, sample, BAM_EXT)


def has_normal(sample):
    return _has_matched_normal_impl(samples_df, sample)


def get_cnv_unit_ids():
    return _get_cnv_unit_ids_impl(samples_df, BAM_EXT)


def get_bam_by_basename(wildcards):
    return BAM_PATH_MAP[wildcards.bam_basename]


def get_unit_tumor_bam_basename(unit):
    return CNV_UNITS[unit]["tumor_bam"]


def get_unit_sample(unit):
    return CNV_UNITS[unit]["sample"]


def get_unit_vcf_path(unit):
    sample = get_unit_sample(unit)
    return _resolve_sample_vcf_impl(samples_df, sample, VCF_FOLDER)


def get_unit_vcf_stats_path(unit):
    sample = get_unit_sample(unit)
    return _resolve_sample_vcf_stats_impl(samples_df, sample, VCF_FOLDER)


def get_unit_vcf_input(wildcards):
    vcf = get_unit_vcf_path(wildcards.unit)
    return [vcf] if vcf else []


def get_unit_vcf_args(wildcards):
    meta = CNV_UNITS[wildcards.unit]
    vcf = get_unit_vcf_path(wildcards.unit)
    return _build_cnvkit_vcf_args_impl(
        vcf=vcf,
        tumor_id=meta["tumor_vcf_id"],
        normal_id=meta["normal_vcf_id"],
        min_variant_depth=CNVKIT_MIN_VARIANT_DEPTH,
    )


def get_access_bed(wildcards=None):
    """Return path to accessible regions BED (pre-supplied, blacklist-subtracted, or generated)."""
    if REF_ACCESS_BED and not REF_BLACKLIST_BED:
        return REF_ACCESS_BED
    return os.path.join(CNVKIT_REF_DIR, "access.bed")


def get_targets_bed_for_binning(wildcards=None):
    """Return path to capture targets BED (blacklist-subtracted when blacklist_bed is set)."""
    if REF_BLACKLIST_BED:
        return os.path.join(CNVKIT_REF_DIR, "targets.clean.bed")
    return TARGETS_BED


def get_cnvkit_reference(wildcards=None):
    """Return pre-built CNVkit reference or pipeline-generated pooled_reference.cnn."""
    if CNVKIT_REFERENCE_CNN:
        return CNVKIT_REFERENCE_CNN
    return os.path.join(CNVKIT_REF_DIR, "pooled_reference.cnn")


def get_purecn_units():
    """Return somatic CNV units eligible for PureCN (must have a VCF and purecn.enabled=true)."""
    if not PURECN_ENABLED:
        return []
    return [u for u in get_cnv_unit_ids() if bool(get_unit_vcf_path(u))]


def get_purecn_mapping_bias_path():
    """Return path to PureCN mapping_bias.rds if configured or generated."""
    if PURECN_MAPPING_BIAS:
        return PURECN_MAPPING_BIAS
    if PURECN_NORMAL_PANEL_VCF:
        return os.path.join(PURECN_REF_DIR, f"mapping_bias_{PURECN_ASSAY}_{PURECN_GENOME}.rds")
    return ""


def get_purecn_normaldb_path():
    """Return path to PureCN normalDB RDS in native mode."""
    if PURECN_NORMALDB:
        return PURECN_NORMALDB
    return os.path.join(PURECN_REF_DIR, f"normalDB_{PURECN_GENOME}.rds")


def get_final_outputs():
    """Return all expected final output files for the cohort."""
    units = get_cnv_unit_ids()
    outputs = []
    if units:
        outputs.extend(expand(os.path.join(CNVKIT_CALLS_DIR, "{unit}.call.cns"), unit=units))
        outputs.extend(expand(os.path.join(CNVKIT_CALLS_DIR, "{unit}.bintest.cns"), unit=units))
        outputs.extend(expand(os.path.join(CNVKIT_EXPORTS_DIR, "{unit}.seg"), unit=units))
        outputs.extend(expand(os.path.join(CNVKIT_EXPORTS_DIR, "{unit}.cnv.bed"), unit=units))
        outputs.extend(expand(os.path.join(CNVKIT_EXPORTS_DIR, "{unit}.cnv.vcf.gz"), unit=units))
        outputs.extend(
            expand(os.path.join(CNVKIT_EXPORTS_DIR, "{unit}.cnv.vcf.gz.tbi"), unit=units)
        )
        outputs.extend(
            expand(os.path.join(QC_DIR, "tumor_metrics", "{unit}.metrics.tsv"), unit=units)
        )
        outputs.extend(
            expand(
                os.path.join(CNVKIT_PLOTS_DIR, "scatter_genome", "{unit}.genome.pdf"),
                unit=units,
            )
        )
        outputs.extend(
            expand(
                os.path.join(CNVKIT_PLOTS_DIR, "diagram", "{unit}.diagram.pdf"),
                unit=units,
            )
        )
        if REF_FLAT and CNVKIT_GENES:
            outputs.extend(
                expand(
                    os.path.join(CNVKIT_PLOTS_DIR, "scatter_gene", "{unit}", "{gene}.pdf"),
                    unit=units,
                    gene=CNVKIT_GENES,
                )
            )
        outputs.append(os.path.join(CNVKIT_PLOTS_DIR, "heatmap_cohort.pdf"))
        outputs.append(os.path.join(QC_DIR, "pipeline_summary.tsv"))
        outputs.append(os.path.join(QC_DIR, "normal_qc_summary.tsv"))
    elif get_normal_bams():
        outputs.append(get_cnvkit_reference())
    return outputs
