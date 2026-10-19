"""Pure helper functions for the sm-cnv-calling workflow.

All functions are free of Snakemake imports so they can be unit-tested directly.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

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


def resolve_baf_vcf(samples_df: pd.DataFrame, sample: str, vcf_folder: str = "") -> str:
    """Resolve BAF VCF path for CNVkit from samples.tsv ('germline_vcf', fallback 'vcf') or paths.vcf_folder."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    for col in ("germline_vcf", "vcf"):
        val = str(row.get(col, "")).strip()
        if val and val != "." and pd.notna(row.get(col)):
            return val
    folder = str(vcf_folder).strip()
    if folder and folder != ".":
        return os.path.join(folder, "mutect2", "filtered", f"{sample}.filtered.vcf.gz")
    return ""


def resolve_purecn_vcf(samples_df: pd.DataFrame, sample: str, vcf_folder: str = "") -> str:
    """Resolve Mutect2 somatic VCF path for PureCN from samples.tsv ('vcf') or paths.vcf_folder.

    PureCN requires a Mutect2 tumor VCF (somatic + germline variants with tumor AD/AF and
    Mutect2 annotations); germline_vcf is intentionally ignored. If no Mutect2 VCF is
    available, logs the reason and returns '' so PureCN is skipped for that sample.
    """
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    val = str(row.get("vcf", "")).strip()
    if val and val != "." and pd.notna(row.get("vcf")):
        return val
    folder = str(vcf_folder).strip()
    if folder and folder != ".":
        return os.path.join(folder, "mutect2", "filtered", f"{sample}.filtered.vcf.gz")
    germline_val = str(row.get("germline_vcf", "")).strip()
    if germline_val and germline_val != "." and pd.notna(row.get("germline_vcf")):
        logger.warning(
            "Skipping PureCN for sample '%s': germline_vcf='%s' is set but no Mutect2 VCF "
            "('vcf' column or 'paths.vcf_folder') was provided.",
            sample,
            germline_val,
        )
    else:
        logger.info(
            "Skipping PureCN for sample '%s': no Mutect2 VCF ('vcf' column or 'paths.vcf_folder') "
            "was provided.",
            sample,
        )
    return ""


def resolve_sample_vcf(samples_df: pd.DataFrame, sample: str, vcf_folder: str = "") -> str:
    """Backward-compatible alias for resolve_baf_vcf()."""
    return resolve_baf_vcf(samples_df, sample, vcf_folder)


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
    """Build CNVkit VCF BAF CLI flags (-v, -i, -n, --min-variant-depth)."""
    vcf_clean = str(vcf).strip()
    if not vcf_clean or vcf_clean == ".":
        return ""
    parts = [f"-v {vcf_clean}"]
    if tumor_id and str(tumor_id).strip() not in ("", "."):
        parts.append(f"-i {str(tumor_id).strip()}")
    if normal_id and str(normal_id).strip() not in ("", "."):
        parts.append(f"-n {str(normal_id).strip()}")
    if min_variant_depth is not None and int(min_variant_depth) > 0:
        parts.append(f"--min-variant-depth {int(min_variant_depth)}")
    return " ".join(parts)


def extract_reference_regions(
    reference_cnn: str | Path,
    targets_bed: str | Path,
    antitargets_bed: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract target and antitarget BED4 intervals from a pre-built CNVkit reference .cnn."""
    df = pd.read_csv(reference_cnn, sep="\t", dtype=str, keep_default_na=False)
    required = {"chromosome", "start", "end", "gene"}
    if not required.issubset(df.columns):
        raise ValueError(f"Reference CNN {reference_cnn} must contain columns {sorted(required)}.")
    cols = ["chromosome", "start", "end", "gene"]
    is_bg = df["gene"].isin(("Antitarget", "Background"))
    targets_df = df.loc[~is_bg, cols].reset_index(drop=True)
    antitargets_df = df.loc[is_bg, cols].reset_index(drop=True)

    t_path = Path(targets_bed)
    a_path = Path(antitargets_bed)
    t_path.parent.mkdir(parents=True, exist_ok=True)
    a_path.parent.mkdir(parents=True, exist_ok=True)
    targets_df.to_csv(t_path, sep="\t", header=False, index=False)
    antitargets_df.to_csv(a_path, sep="\t", header=False, index=False)
    return targets_df, antitargets_df


def identify_clean_normals(
    metrics_file: str | Path,
    coverage_dir: str | Path,
    output_list: str | Path,
    threshold: float = 0.3,
) -> list[str]:
    """Filter normal QC metrics by bivar <= threshold and resolve verified .cnn paths.

    Writes both .targetcoverage.cnn and .antitargetcoverage.cnn paths for each
    passing normal sample to output_list and returns the list of paths.
    """
    metrics_path = Path(metrics_file)
    cov_dir = Path(coverage_dir)
    out_path = Path(output_list)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not metrics_path.exists() or metrics_path.stat().st_size == 0:
        out_path.write_text("", encoding="utf-8")
        return []

    df = pd.read_csv(metrics_path, sep="\t")
    if "sample" not in df.columns or "bivar" not in df.columns:
        raise ValueError(f"Metrics file {metrics_file} must contain 'sample' and 'bivar' columns.")

    df["bivar"] = pd.to_numeric(df["bivar"], errors="coerce")
    clean_df = df[df["bivar"].notna() & (df["bivar"] <= float(threshold))]

    resolved_paths: list[str] = []
    seen_basenames: set[str] = set()
    for raw_sample in clean_df["sample"].astype(str).tolist():
        bn = os.path.basename(raw_sample.strip())
        for suffix in (".targetcoverage.cnn", ".antitargetcoverage.cnn", ".cnn"):
            if bn.endswith(suffix):
                bn = bn[: -len(suffix)]
                break
        if bn in seen_basenames:
            continue
        seen_basenames.add(bn)

        target_cnn = cov_dir / f"{bn}.targetcoverage.cnn"
        antitarget_cnn = cov_dir / f"{bn}.antitargetcoverage.cnn"
        if not target_cnn.is_file():
            raise FileNotFoundError(f"Clean normal target coverage file not found: {target_cnn}")
        if not antitarget_cnn.is_file():
            raise FileNotFoundError(
                f"Clean normal antitarget coverage file not found: {antitarget_cnn}"
            )
        resolved_paths.append(str(target_cnn))
        resolved_paths.append(str(antitarget_cnn))

    content = "\n".join(resolved_paths) + ("\n" if resolved_paths else "")
    out_path.write_text(content, encoding="utf-8")
    return resolved_paths


def parse_cnvkit_sex(sex_file: str | Path = "", sample_sex_override: str = "") -> str:
    """Resolve chromosomal sex ('male' or 'female') from override or cnvkit sex output."""
    override = normalize_sex(sample_sex_override)
    if override:
        return override

    if sex_file and Path(sex_file).is_file() and Path(sex_file).stat().st_size > 0:
        try:
            df = pd.read_csv(sex_file, sep="\t")
            if not df.empty and "sex" in df.columns:
                inferred = normalize_sex(str(df.iloc[0]["sex"]))
                if inferred:
                    return inferred
        except Exception:
            pass
    return "female"


def consolidate_purity_ploidy(
    unit: str,
    sample: str,
    output_file: str | Path,
    purecn_csv: str | Path = "",
    purecn_status_file: str | Path = "",
    sample_purity: float | str | None = None,
    sample_ploidy: float | str | None = None,
    purecn_attempted: bool = False,
) -> dict[str, Any]:
    """Consolidate tumor purity and ploidy from PureCN or explicit sample metadata.

    Never fabricates a placeholder purity: when neither PureCN nor explicit sample
    metadata provides a valid purity, records source='unestimated' so cnvkit call
    falls back to threshold mode.
    """
    purity_val: float | None = None
    ploidy_val: float | None = None
    source = "unestimated"
    status = "attempted" if purecn_attempted else "skipped"

    if purecn_status_file and Path(purecn_status_file).is_file():
        raw_status = Path(purecn_status_file).read_text(encoding="utf-8").strip().lower()
        if raw_status in ("succeeded", "failed", "skipped"):
            status = raw_status

    if purecn_csv and Path(purecn_csv).is_file() and Path(purecn_csv).stat().st_size > 0:
        try:
            df = pd.read_csv(purecn_csv)
            if not df.empty and "Purity" in df.columns:
                p_top = float(df.iloc[0]["Purity"])
                pl_top = (
                    float(df.iloc[0]["Ploidy"])
                    if "Ploidy" in df.columns and pd.notna(df.iloc[0]["Ploidy"])
                    else 2.0
                )
                if pd.notna(p_top) and 0.0 < p_top <= 1.0 and pl_top > 0.0:
                    purity_val = round(p_top, 4)
                    ploidy_val = round(pl_top, 4)
                    source = "PureCN"
                    status = "succeeded"
                else:
                    status = "failed"
            else:
                status = "failed"
        except Exception:
            status = "failed"
    elif purecn_attempted:
        status = "failed"

    if purity_val is None and sample_purity is not None:
        raw_p = str(sample_purity).strip()
        if raw_p and raw_p not in ("", ".", "None", "nan"):
            p_float = float(raw_p)
            if not (0.0 < p_float <= 1.0):
                raise ValueError(f"Invalid explicit sample purity: {p_float}")
            pl_float = 2.0
            if sample_ploidy is not None:
                raw_pl = str(sample_ploidy).strip()
                if raw_pl and raw_pl not in ("", ".", "None", "nan"):
                    pl_float = float(raw_pl)
                    if pl_float <= 0.0:
                        raise ValueError(f"Invalid explicit sample ploidy: {pl_float}")
            purity_val = round(p_float, 4)
            ploidy_val = round(pl_float, 4)
            source = "sample_sheet"

    record = {
        "unit": unit,
        "sample": sample,
        "purity": "" if purity_val is None else purity_val,
        "ploidy": "" if ploidy_val is None else ploidy_val,
        "source": source,
        "status": status,
    }
    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([record]).to_csv(out_path, sep="\t", index=False)
    return record


def build_cnvkit_call_args(
    purity_ploidy_tsv: str | Path,
    sex: str = "",
    male_reference: bool = False,
    call_method: str = "clonal",
    call_filter: str = "ci",
    drop_low_coverage: bool = True,
    vcf_args: str = "",
    extra: str = "",
) -> str:
    """Build CLI flags for cnvkit.py call from consolidated purity/ploidy and sex."""
    df = pd.read_csv(purity_ploidy_tsv, sep="\t", dtype=str, keep_default_na=False)
    row = df.iloc[0] if not df.empty else {}
    purity_str = str(row.get("purity", "")).strip()
    ploidy_str = str(row.get("ploidy", "")).strip()
    source = str(row.get("source", "unestimated")).strip()

    args: list[str] = []
    if call_filter and call_filter != "none":
        args.append(f"--filter {call_filter}")
    if drop_low_coverage:
        args.append("--drop-low-coverage")

    if (
        call_method == "clonal"
        and source in ("PureCN", "sample_sheet")
        and purity_str not in ("", ".")
    ):
        args.append("-m clonal")
        args.append(f"--purity {purity_str}")
        if ploidy_str not in ("", "."):
            args.append(f"--ploidy {round(float(ploidy_str))}")
    elif call_method == "none":
        args.append("-m none")
    else:
        args.append("-m threshold")

    if male_reference:
        args.append("-y")
    norm_sex = normalize_sex(sex)
    if norm_sex:
        args.append(f"-x {norm_sex}")
    if vcf_args:
        args.append(vcf_args.strip())
    if extra and extra.strip():
        args.append(extra.strip())
    return " ".join(args)


def build_pipeline_summary(
    units_meta: dict[str, dict[str, Any]],
    purity_files: list[str],
    sex_files: list[str],
    tumor_metrics_files: list[str],
    normal_metrics_file: str | None,
    clean_normals_file: str | None,
    output_summary_tsv: str | Path,
    output_normal_qc_tsv: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create cohort pipeline_summary.tsv and normal_qc_summary.tsv reports."""
    purity_map: dict[str, dict[str, str]] = {}
    for pf in purity_files:
        if Path(pf).is_file() and Path(pf).stat().st_size > 0:
            pdf = pd.read_csv(pf, sep="\t", dtype=str, keep_default_na=False)
            for _, r in pdf.iterrows():
                purity_map[str(r["unit"])] = {
                    "purity": str(r.get("purity", "")),
                    "ploidy": str(r.get("ploidy", "")),
                    "purity_source": str(r.get("source", "unestimated")),
                    "purecn_status": str(r.get("status", "skipped")),
                }

    sex_map: dict[str, str] = {}
    for sf in sex_files:
        p = Path(sf)
        if p.is_file():
            unit_key = p.name.replace(".sex.txt", "")
            sex_map[unit_key] = p.read_text(encoding="utf-8").strip()

    metrics_map: dict[str, dict[str, Any]] = {}
    for mf in tumor_metrics_files:
        p = Path(mf)
        if p.is_file() and p.stat().st_size > 0:
            unit_key = p.name.replace(".metrics.tsv", "")
            try:
                mdf = pd.read_csv(p, sep="\t")
                if not mdf.empty:
                    r0 = mdf.iloc[0]
                    metrics_map[unit_key] = {
                        "segments": r0.get("segments", ""),
                        "stdev": r0.get("stdev", ""),
                        "mad": r0.get("mad", ""),
                        "iqr": r0.get("iqr", ""),
                        "bivar": r0.get("bivar", ""),
                    }
            except Exception:
                pass

    summary_rows: list[dict[str, Any]] = []
    for unit_id, meta in units_meta.items():
        pinfo = purity_map.get(unit_id, {})
        minfo = metrics_map.get(unit_id, {})
        summary_rows.append(
            {
                "unit": unit_id,
                "sample": meta["sample"],
                "analysis_type": meta["analysis_type"],
                "tumor_bam": meta["tumor_bam"],
                "normal_bam": meta["normal_bam"],
                "sex": sex_map.get(unit_id, meta.get("sex", "")),
                "purity": pinfo.get("purity", ""),
                "ploidy": pinfo.get("ploidy", ""),
                "purity_source": pinfo.get("purity_source", "unestimated"),
                "purecn_status": pinfo.get("purecn_status", "skipped"),
                "segments": minfo.get("segments", ""),
                "stdev": minfo.get("stdev", ""),
                "mad": minfo.get("mad", ""),
                "iqr": minfo.get("iqr", ""),
                "bivar": minfo.get("bivar", ""),
            }
        )

    summary_cols = [
        "unit",
        "sample",
        "analysis_type",
        "tumor_bam",
        "normal_bam",
        "sex",
        "purity",
        "ploidy",
        "purity_source",
        "purecn_status",
        "segments",
        "stdev",
        "mad",
        "iqr",
        "bivar",
    ]
    summary_df = pd.DataFrame(summary_rows, columns=summary_cols)
    out_sum = Path(output_summary_tsv)
    out_sum.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(out_sum, sep="\t", index=False)

    clean_set: set[str] = set()
    if clean_normals_file and Path(clean_normals_file).is_file():
        for line in Path(clean_normals_file).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line.endswith(".targetcoverage.cnn"):
                clean_set.add(os.path.basename(line)[: -len(".targetcoverage.cnn")])

    normal_rows: list[dict[str, Any]] = []
    if (
        normal_metrics_file
        and Path(normal_metrics_file).is_file()
        and Path(normal_metrics_file).stat().st_size > 0
    ):
        try:
            ndf = pd.read_csv(normal_metrics_file, sep="\t")
            for _, nrow in ndf.iterrows():
                raw_s = os.path.basename(str(nrow.get("sample", "")).strip())
                for sfx in (".targetcoverage.cnn", ".antitargetcoverage.cnn", ".cnn"):
                    if raw_s.endswith(sfx):
                        raw_s = raw_s[: -len(sfx)]
                        break
                normal_rows.append(
                    {
                        "normal_bam": raw_s,
                        "stdev": nrow.get("stdev", ""),
                        "mad": nrow.get("mad", ""),
                        "iqr": nrow.get("iqr", ""),
                        "bivar": nrow.get("bivar", ""),
                        "included_in_pon": raw_s in clean_set,
                    }
                )
        except Exception:
            pass

    normal_cols = ["normal_bam", "stdev", "mad", "iqr", "bivar", "included_in_pon"]
    normal_df = pd.DataFrame(normal_rows, columns=normal_cols)
    out_norm = Path(output_normal_qc_tsv)
    out_norm.parent.mkdir(parents=True, exist_ok=True)
    normal_df.to_csv(out_norm, sep="\t", index=False)

    return summary_df, normal_df
