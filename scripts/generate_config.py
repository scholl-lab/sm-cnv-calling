#!/usr/bin/env python3
"""Generate config/samples.tsv and config/config.yaml for sm-cnv-calling.

Scans a BAM directory (e.g. from sm-alignment) or imports an sm-calling samples.tsv,
classifies tumor/normal roles, pairs tumors with normals, and writes config files.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import pandas as pd

DEFAULT_BAM_EXTENSION = ".merged.dedup.bqsr.bam"

TUMOR_PATTERNS = [
    re.compile(r"[._-](?:tumor|tumour|FFPE|met|metastatic|primary|tum)$", re.IGNORECASE),
    re.compile(r"[._-]T[A-Z]?\d*$", re.IGNORECASE),
    re.compile(r"[._-]L\d*$", re.IGNORECASE),
]

NORMAL_PATTERNS = [
    re.compile(r"[._-](?:normal|blood|buffy|germline|pbmc|ctrl)$", re.IGNORECASE),
    re.compile(r"[._-]N\d*$", re.IGNORECASE),
]

ROLE_SUFFIX_RE = re.compile(
    r"[._-](?:tumor|tumour|FFPE|met|metastatic|primary|tum|"
    r"normal|blood|buffy|germline|pbmc|ctrl|"
    r"T[A-Z]?\d*|N\d*|L\d*)$",
    re.IGNORECASE,
)


def guess_role(basename: str) -> str:
    """Classify a BAM basename as 'tumor', 'normal', or 'unknown'."""
    for pat in TUMOR_PATTERNS:
        if pat.search(basename):
            return "tumor"
    for pat in NORMAL_PATTERNS:
        if pat.search(basename):
            return "normal"
    return "unknown"


def extract_patient_id(basename: str) -> str:
    """Extract patient identifier by stripping tumor/normal role suffix."""
    return ROLE_SUFFIX_RE.sub("", basename)


def discover_bam_files(bam_folder: str | Path, bam_extension: str) -> list[dict[str, Any]]:
    """Scan a directory for BAM files matching bam_extension."""
    bam_dir = Path(bam_folder)
    if not bam_dir.is_dir():
        return []
    entries: list[dict[str, Any]] = []
    for path in sorted(bam_dir.iterdir()):
        if path.is_file() and path.name.endswith(bam_extension):
            bn = path.name[: -len(bam_extension)]
            entries.append(
                {
                    "basename": bn,
                    "path": str(path),
                    "role": guess_role(bn),
                }
            )
    return entries


def build_sample_rows(
    bam_entries: list[dict[str, Any]],
    analysis_mode: str = "auto",
) -> list[dict[str, str]]:
    """Pair tumor and normal BAM entries into sm-cnv-calling sample sheet rows."""
    groups: dict[str, dict[str, list[str]]] = {}
    for entry in bam_entries:
        pid = extract_patient_id(entry["basename"])
        groups.setdefault(pid, {"tumors": [], "normals": []})
        if entry["role"] == "normal":
            groups[pid]["normals"].append(entry["basename"])
        else:
            groups[pid]["tumors"].append(entry["basename"])

    rows: list[dict[str, str]] = []
    for pid, grp in groups.items():
        tumors = grp["tumors"]
        normals = grp["normals"]
        matched_normal = normals[0] if normals else "."

        if tumors:
            if analysis_mode in ("auto", "tumor_normal") and matched_normal != ".":
                rows.append(
                    {
                        "sample": f"{pid}_TN",
                        "tumor_bam": ";".join(tumors),
                        "normal_bam": matched_normal,
                        "analysis_type": "tumor_normal",
                        "in_pon": "yes",
                    }
                )
            elif analysis_mode in ("auto", "tumor_only"):
                for idx, t_bn in enumerate(tumors):
                    sfx = f"_{idx + 1}" if len(tumors) > 1 else ""
                    rows.append(
                        {
                            "sample": f"{pid}_To{sfx}",
                            "tumor_bam": t_bn,
                            "normal_bam": ".",
                            "analysis_type": "tumor_only",
                            "in_pon": "no",
                        }
                    )
            for extra_n in normals[1:] if matched_normal != "." else normals:
                rows.append(
                    {
                        "sample": f"{extra_n}_G",
                        "tumor_bam": extra_n,
                        "normal_bam": ".",
                        "analysis_type": "germline",
                        "in_pon": "yes",
                    }
                )
        else:
            for n_bn in normals:
                rows.append(
                    {
                        "sample": f"{pid}_G" if len(normals) == 1 else f"{n_bn}_G",
                        "tumor_bam": n_bn,
                        "normal_bam": ".",
                        "analysis_type": "germline",
                        "in_pon": "yes",
                    }
                )
    return rows


def render_config_yaml(
    samples_path: str = "config/samples.tsv",
    bam_folder: str = "results/exomes/bqsr",
    vcf_folder: str = "results/exomes/variant_calls",
    output_folder: str = "results/exomes/cnv_calls",
    bam_extension: str = DEFAULT_BAM_EXTENSION,
    genome_fasta: str = "resources/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    targets_bed: str = "resources/targets/capture_targets.bed",
) -> str:
    """Render a schema-valid config.yaml string with GRCh38 defaults."""
    return f"""# sm-cnv-calling configuration (generated by scripts/generate_config.py)

ref:
  genome: "{genome_fasta}"
  build: "GRCh38"
  targets_bed: "{targets_bed}"
  baits_bed: ""
  access_bed: ""
  blacklist_bed: ""
  ref_flat: ""
  mappability_bw: ""
  reptiming_bw: ""

paths:
  samples: "{samples_path}"
  bam_folder: "{bam_folder}"
  vcf_folder: "{vcf_folder}"
  output_folder: "{output_folder}"
  log_subdir: "logs"

bam:
  file_extension: "{bam_extension}"

cnvkit:
  reference_cnn: ""
  autobin: false
  target_avg_size: 200
  antitarget_avg_size: 150000
  antitarget_min_size: 500
  access_min_gap_size: 5000
  pon_bivar_threshold: 0.3
  male_reference: false
  segment_method: "cbs"
  segment_threshold: 1.0e-4
  drop_low_coverage: true
  call_method: "clonal"
  call_filter: "ci"
  min_variant_depth: 20
  genes_of_interest:
    - "TP53"
    - "MYC"
    - "EGFR"
    - "ERBB2"

purecn:
  enabled: true
  mode: "cnvkit"
  genome: "hg38"
  assay: "exome"
  fun_segmentation: "none"
  normal_panel_vcf: ""
  normaldb: ""
  mapping_bias: ""
  snp_blacklist: ""
  seed: 123
  postoptimize: true
  extra: ""

params:
  cnvkit:
    segment_method: "cbs"
    segment_threshold: 1.0e-4
    access_extra: ""
    target_extra: ""
    antitarget_extra: ""
    autobin_extra: ""
    coverage_extra: ""
    reference_extra: ""
    fix_extra: ""
    segment_extra: ""
    segmetrics_extra: ""
    call_extra: ""
    bintest_extra: ""
    scatter_extra: ""
    diagram_extra: ""
    heatmap_extra: ""
  purecn:
    intervalfile_extra: ""
    coverage_extra: ""
    normaldb_extra: ""
    run_extra: ""
"""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate config/samples.tsv and config/config.yaml for sm-cnv-calling."
    )
    parser.add_argument("--bam-folder", required=True, help="Directory with input BAM files.")
    parser.add_argument(
        "--bam-extension",
        default=DEFAULT_BAM_EXTENSION,
        help=f"BAM filename extension (default: {DEFAULT_BAM_EXTENSION}).",
    )
    parser.add_argument(
        "--analysis-mode",
        choices=["auto", "tumor_normal", "tumor_only"],
        default="auto",
        help="Pairing mode for discovered BAMs.",
    )
    parser.add_argument("--output", default="config/samples.tsv", help="Output samples.tsv path.")
    parser.add_argument(
        "--config-template",
        action="store_true",
        help="Also write config.yaml to --config-output.",
    )
    parser.add_argument(
        "--config-output", default="config/config.yaml", help="Output config.yaml path."
    )
    parser.add_argument(
        "--vcf-folder",
        default="results/exomes/variant_calls",
        help="Path to sm-calling output folder.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/exomes/cnv_calls",
        help="Pipeline output folder.",
    )
    parser.add_argument(
        "--targets-bed",
        default="resources/targets/capture_targets.bed",
        help="Capture targets BED path.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print without writing files.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing output files.")

    args = parser.parse_args(argv)
    entries = discover_bam_files(args.bam_folder, args.bam_extension)
    rows = build_sample_rows(entries, args.analysis_mode)
    df = pd.DataFrame(
        rows, columns=["sample", "tumor_bam", "normal_bam", "analysis_type", "in_pon"]
    )

    if args.dry_run:
        print(df.to_csv(sep="\t", index=False))
        return 0

    out_samples = Path(args.output)
    if out_samples.exists() and not args.force:
        raise FileExistsError(f"{out_samples} already exists (use --force to overwrite).")
    out_samples.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_samples, sep="\t", index=False)

    if args.config_template:
        out_cfg = Path(args.config_output)
        if out_cfg.exists() and not args.force:
            raise FileExistsError(f"{out_cfg} already exists (use --force to overwrite).")
        out_cfg.parent.mkdir(parents=True, exist_ok=True)
        out_cfg.write_text(
            render_config_yaml(
                samples_path=str(out_samples),
                bam_folder=args.bam_folder,
                vcf_folder=args.vcf_folder,
                output_folder=args.output_dir,
                bam_extension=args.bam_extension,
                targets_bed=args.targets_bed,
            ),
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
