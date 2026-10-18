"""Pure helper functions for the sm-cnv-calling workflow.

All functions are free of Snakemake imports so they can be unit-tested directly.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import pandas as pd

DEFAULT_BAM_EXT = ".merged.dedup.bqsr.bam"

ANALYSIS_TYPE_ALIASES: dict[str, str] = {
    "tumor_normal": "tumor_normal",
    "tumor_only": "tumor_only",
    "germline": "germline",
    "TvsN": "tumor_normal",
    "To": "tumor_only",
}

SEX_ALIASES: dict[str, str] = {
    "m": "male",
    "male": "male",
    "y": "male",
    "f": "female",
    "female": "female",
    "x": "female",
}


def normalize_analysis_type(analysis_type: str) -> str:
    """Normalize legacy analysis_type aliases ('TvsN'/'To') to canonical sm-calling names."""
    val = str(analysis_type).strip()
    if val not in ANALYSIS_TYPE_ALIASES:
        raise ValueError(
            f"Invalid analysis_type '{analysis_type}'. "
            f"Expected one of {sorted(ANALYSIS_TYPE_ALIASES.keys())}."
        )
    return ANALYSIS_TYPE_ALIASES[val]


def normalize_sex(sex: str | None) -> str:
    """Normalize chromosomal sex string to 'male', 'female', or ''."""
    if sex is None or pd.isna(sex):
        return ""
    val = str(sex).strip().lower()
    if val in ("", "."):
        return ""
    if val not in SEX_ALIASES:
        raise ValueError(
            f"Invalid sex value '{sex}'. Expected 'male'/'female', 'M'/'F', 'x'/'y', or '.'."
        )
    return SEX_ALIASES[val]


def strip_bam_extension(bam_entry: str, bam_ext: str = DEFAULT_BAM_EXT) -> str:
    """Return clean BAM basename by stripping directory and bam_ext (or .bam)."""
    name = os.path.basename(str(bam_entry).strip())
    if bam_ext and name.endswith(bam_ext):
        return name[: -len(bam_ext)]
    if name.endswith(".bam"):
        return name[:-4]
    return name


def resolve_bam_path(bam_entry: str, bam_folder: str, bam_ext: str = DEFAULT_BAM_EXT) -> str:
    """Resolve a BAM basename or file path to a full path."""
    raw = str(bam_entry).strip()
    if not raw or raw == ".":
        raise ValueError("Cannot resolve empty or '.' BAM entry.")
    if os.path.isabs(raw) or "/" in raw or "\\" in raw:
        if raw.endswith(".bam"):
            return raw
        return raw + bam_ext
    clean = strip_bam_extension(raw, bam_ext)
    return os.path.join(bam_folder, clean + bam_ext)


def validate_samples_df(samples_df: pd.DataFrame) -> None:
    """Validate sample sheet invariants beyond per-row JSON schema checks."""
    if "sample" not in samples_df.columns:
        raise ValueError("Sample sheet is missing required 'sample' column.")
    if samples_df["sample"].duplicated().any() or samples_df.index.duplicated().any():
        dups = sorted(
            samples_df.loc[samples_df["sample"].duplicated(), "sample"].astype(str).unique()
        )
        raise ValueError(f"Duplicate sample IDs in samples.tsv: {dups}")

    for _, row in samples_df.iterrows():
        sample = row.get("sample", "")
        tumor = row.get("tumor_bam", "")
        normal = row.get("normal_bam", "")
        raw_type = row.get("analysis_type", "")

        if not isinstance(sample, str) or not isinstance(tumor, str):
            raise ValueError(
                "Sample and BAM identifiers must be strings (read samples.tsv with dtype=str)"
            )
        sample_str = sample.strip()
        tumor_parts = [b.strip() for b in tumor.split(";")]
        if (
            not sample_str
            or sample_str == "."
            or not tumor.strip()
            or any(not b or b == "." for b in tumor_parts)
        ):
            raise ValueError(f"Sample '{sample}' has empty or '.' sample/tumor_bam")

        if len(set(tumor_parts)) != len(tumor_parts):
            raise ValueError(
                f"Sample '{sample}' has duplicate BAM entries in semicolon-delimited tumor_bam: {tumor}"
            )

        analysis_type = normalize_analysis_type(str(raw_type))

        has_normal_bam = bool(normal) and normal != "." and pd.notna(normal)
        if has_normal_bam and not isinstance(normal, str):
            raise ValueError("normal_bam must be a string (read samples.tsv with dtype=str)")
        if has_normal_bam and ";" in str(normal):
            raise ValueError(
                f"Sample '{sample}' has multiple semicolon-delimited entries in normal_bam: {normal}"
            )
        if analysis_type == "tumor_normal" and not has_normal_bam:
            raise ValueError(
                f"Sample '{sample}' has analysis_type='{raw_type}' but missing normal_bam"
            )
        if analysis_type != "tumor_normal" and has_normal_bam:
            raise ValueError(
                f"Sample '{sample}' has analysis_type='{raw_type}' "
                f"but normal_bam='{normal}' is provided"
            )

        raw_in_pon = str(row.get("in_pon", "")).strip().lower()
        if raw_in_pon not in ("", ".", "yes", "no", "true", "false", "1", "0"):
            raise ValueError(
                f"Sample '{sample}' has invalid in_pon='{row.get('in_pon')}' "
                "(expected 'yes'/'no' or '.')"
            )

        normalize_sex(row.get("sex", ""))

        raw_purity = str(row.get("purity", "")).strip()
        if raw_purity and raw_purity != ".":
            try:
                pval = float(raw_purity)
            except ValueError as exc:
                raise ValueError(
                    f"Sample '{sample}' has non-numeric purity='{raw_purity}'"
                ) from exc
            if not (0.0 < pval <= 1.0):
                raise ValueError(f"Sample '{sample}' has purity={pval} outside valid range (0, 1]")

        raw_ploidy = str(row.get("ploidy", "")).strip()
        if raw_ploidy and raw_ploidy != ".":
            try:
                plval = float(raw_ploidy)
            except ValueError as exc:
                raise ValueError(
                    f"Sample '{sample}' has non-numeric ploidy='{raw_ploidy}'"
                ) from exc
            if plval <= 0.0:
                raise ValueError(f"Sample '{sample}' has ploidy={plval} (must be > 0)")


def get_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return sorted list of unique sample identifiers."""
    return sorted(samples_df["sample"].unique().tolist())


def get_tumor_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return sorted list of somatic sample IDs (tumor_only or tumor_normal)."""
    somatic: list[str] = []
    for _, row in samples_df.iterrows():
        atype = normalize_analysis_type(str(row.get("analysis_type", "")))
        if atype in ("tumor_only", "tumor_normal"):
            somatic.append(str(row["sample"]))
    return sorted(set(somatic))


def get_germline_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return sorted list of germline sample IDs."""
    germline: list[str] = []
    for _, row in samples_df.iterrows():
        atype = normalize_analysis_type(str(row.get("analysis_type", "")))
        if atype == "germline":
            germline.append(str(row["sample"]))
    return sorted(set(germline))


def get_tumor_bam_basenames(
    samples_df: pd.DataFrame, sample: str, bam_ext: str = DEFAULT_BAM_EXT
) -> list[str]:
    """Return tumor BAM basenames for a sample (splitting ';'-delimited multi-region entries)."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    raw = str(row.get("tumor_bam", ""))
    return [
        strip_bam_extension(b.strip(), bam_ext)
        for b in raw.split(";")
        if b.strip() and b.strip() != "."
    ]


def get_normal_bam_basenames(
    samples_df: pd.DataFrame, bam_ext: str = DEFAULT_BAM_EXT
) -> list[str]:
    """Return sorted basenames of normal BAMs for PoN and PureCN NormalDB."""
    has_in_pon = "in_pon" in samples_df.columns
    normals: set[str] = set()
    for _, row in samples_df.iterrows():
        raw_in_pon = str(row.get("in_pon", "")).strip().lower()
        explicitly_excluded = has_in_pon and raw_in_pon in ("no", "false", "0")
        explicitly_included = raw_in_pon in ("yes", "true", "1")

        atype = normalize_analysis_type(str(row.get("analysis_type", "")))
        normal = row.get("normal_bam", "")
        if atype == "tumor_normal" and normal and normal != "." and pd.notna(normal):
            if not explicitly_excluded:
                normals.add(strip_bam_extension(str(normal), bam_ext))
        elif atype == "germline" and has_in_pon and explicitly_included:
            for b in str(row.get("tumor_bam", "")).split(";"):
                b = b.strip()
                if b and b != ".":
                    normals.add(strip_bam_extension(b, bam_ext))
    return sorted(normals)


def get_unique_bam_basenames(
    samples_df: pd.DataFrame, bam_ext: str = DEFAULT_BAM_EXT
) -> list[str]:
    """Return all unique BAM basenames required across somatic tumors and PoN normals."""
    bams: set[str] = set(get_normal_bam_basenames(samples_df, bam_ext))
    for _, row in samples_df.iterrows():
        atype = normalize_analysis_type(str(row.get("analysis_type", "")))
        if atype in ("tumor_normal", "tumor_only"):
            for b in str(row.get("tumor_bam", "")).split(";"):
                b = b.strip()
                if b and b != ".":
                    bams.add(strip_bam_extension(b, bam_ext))
        if (
            atype == "tumor_normal"
            and pd.notna(row.get("normal_bam"))
            and str(row.get("normal_bam", "")) not in ("", ".")
        ):
            bams.add(strip_bam_extension(str(row["normal_bam"]), bam_ext))
    return sorted(bams)


def get_bam_path_map(
    samples_df: pd.DataFrame, bam_folder: str, bam_ext: str = DEFAULT_BAM_EXT
) -> dict[str, str]:
    """Build a mapping from clean BAM basename to resolved BAM file path."""
    path_map: dict[str, str] = {}
    for _, row in samples_df.iterrows():
        for raw_b in str(row.get("tumor_bam", "")).split(";"):
            raw_b = raw_b.strip()
            if raw_b and raw_b != ".":
                bn = strip_bam_extension(raw_b, bam_ext)
                path_map[bn] = resolve_bam_path(raw_b, bam_folder, bam_ext)
        normal = row.get("normal_bam", "")
        if pd.notna(normal) and str(normal).strip() not in ("", "."):
            raw_n = str(normal).strip()
            bn = strip_bam_extension(raw_n, bam_ext)
            path_map[bn] = resolve_bam_path(raw_n, bam_folder, bam_ext)
    return path_map


def has_matched_normal(samples_df: pd.DataFrame, sample: str) -> bool:
    """Check whether a sample has a matched normal BAM."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    if normalize_analysis_type(str(row.get("analysis_type", ""))) != "tumor_normal":
        return False
    normal = row.get("normal_bam", "")
    return bool(normal) and normal != "." and pd.notna(normal)


def get_cnv_units(
    samples_df: pd.DataFrame, bam_ext: str = DEFAULT_BAM_EXT
) -> dict[str, dict[str, Any]]:
    """Expand somatic rows into per-tumor-BAM CNV calling units.

    Single-BAM rows use unit ID '{sample}'.
    Multi-region ';'-delimited rows use unit IDs '{sample}__{bam_basename}'.
    """
    units: dict[str, dict[str, Any]] = {}
    for _, row in samples_df.iterrows():
        sample = str(row["sample"]).strip()
        atype = normalize_analysis_type(str(row["analysis_type"]))
        if atype not in ("tumor_normal", "tumor_only"):
            continue

        raw_parts = [
            b.strip() for b in str(row["tumor_bam"]).split(";") if b.strip() and b.strip() != "."
        ]
        basenames = [strip_bam_extension(p, bam_ext) for p in raw_parts]

        normal_raw = str(row.get("normal_bam", "")).strip()
        normal_bn = (
            strip_bam_extension(normal_raw, bam_ext)
            if atype == "tumor_normal" and normal_raw not in ("", ".")
            else ""
        )

        raw_tumor_vcf_ids = [
            x.strip()
            for x in str(row.get("tumor_sample_id_vcf", "")).split(";")
            if x.strip() and x.strip() != "."
        ]
        raw_normal_vcf_id = str(row.get("normal_sample_id_vcf", "")).strip()
        normal_vcf_id = raw_normal_vcf_id if raw_normal_vcf_id not in ("", ".") else normal_bn

        sex = normalize_sex(row.get("sex", ""))

        raw_purity = str(row.get("purity", "")).strip()
        purity: float | None = float(raw_purity) if raw_purity and raw_purity != "." else None
        raw_ploidy = str(row.get("ploidy", "")).strip()
        ploidy: float | None = float(raw_ploidy) if raw_ploidy and raw_ploidy != "." else None

        is_multi = len(basenames) > 1
        for idx, (raw_part, bn) in enumerate(zip(raw_parts, basenames, strict=True)):
            unit_id = f"{sample}__{bn}" if is_multi else sample
            if len(raw_tumor_vcf_ids) == len(basenames):
                tumor_vcf_id = raw_tumor_vcf_ids[idx]
            elif len(raw_tumor_vcf_ids) == 1:
                tumor_vcf_id = raw_tumor_vcf_ids[0]
            else:
                tumor_vcf_id = bn

            units[unit_id] = {
                "unit": unit_id,
                "sample": sample,
                "analysis_type": atype,
                "tumor_bam": bn,
                "tumor_bam_raw": raw_part,
                "normal_bam": normal_bn,
                "normal_bam_raw": normal_raw if normal_bn else "",
                "tumor_vcf_id": tumor_vcf_id,
                "normal_vcf_id": normal_vcf_id if atype == "tumor_normal" else "",
                "sex": sex,
                "purity": purity,
                "ploidy": ploidy,
            }
    return units


def get_cnv_unit_ids(samples_df: pd.DataFrame, bam_ext: str = DEFAULT_BAM_EXT) -> list[str]:
    """Return ordered list of somatic CNV unit IDs."""
    return list(get_cnv_units(samples_df, bam_ext).keys())


def resolve_sample_vcf(samples_df: pd.DataFrame, sample: str, vcf_folder: str = "") -> str:
    """Resolve VCF path for a sample from samples.tsv or paths.vcf_folder."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    for col in ("vcf", "germline_vcf"):
        val = str(row.get(col, "")).strip()
        if val and val != "." and pd.notna(row.get(col)):
            return val
    folder = str(vcf_folder).strip()
    if folder and folder != ".":
        return os.path.join(folder, "mutect2", "filtered", f"{sample}.filtered.vcf.gz")
    return ""


def resolve_sample_vcf_stats(samples_df: pd.DataFrame, sample: str, vcf_folder: str = "") -> str:
    """Resolve Mutect2 .vcf.gz.stats path for a sample if available."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    explicit_stats = str(row.get("vcf_stats", "")).strip()
    if explicit_stats and explicit_stats != "." and pd.notna(row.get("vcf_stats")):
        return explicit_stats
    explicit_vcf = str(row.get("vcf", "")).strip()
    folder = str(vcf_folder).strip()
    if (not explicit_vcf or explicit_vcf == ".") and folder and folder != ".":
        return os.path.join(folder, "mutect2", "merged", f"{sample}.vcf.gz.stats")
    return ""


def build_cnvkit_vcf_args(
    vcf: str = "",
    tumor_id: str = "",
    normal_id: str = "",
    min_variant_depth: int | None = None,
) -> str:
    """Build CNVkit VCF BAF CLI flags (-v, -i, -n, -m)."""
    vcf_clean = str(vcf).strip()
    if not vcf_clean or vcf_clean == ".":
        return ""
    parts = [f"-v {vcf_clean}"]
    if tumor_id and str(tumor_id).strip() not in ("", "."):
        parts.append(f"-i {str(tumor_id).strip()}")
    if normal_id and str(normal_id).strip() not in ("", "."):
        parts.append(f"-n {str(normal_id).strip()}")
    if min_variant_depth is not None and int(min_variant_depth) > 0:
        parts.append(f"-m {int(min_variant_depth)}")
    return " ".join(parts)


