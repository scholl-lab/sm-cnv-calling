"""Shared fixtures and tool resolution for the sm-cnv-calling test suite."""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _augment_environment_path() -> None:
    """Ensure local conda/venv tool directories are on PATH if present."""
    candidate_dirs = [
        Path("/tmp/sm-venv/bin"),
        Path.home() / "miniforge3" / "envs" / "varlens-annotate" / "bin",
        Path.home() / "miniforge3" / "bin",
    ]
    current_parts = os.environ.get("PATH", "").split(os.pathsep)
    prepend = [str(d) for d in candidate_dirs if d.is_dir() and str(d) not in current_parts]
    if prepend:
        os.environ["PATH"] = os.pathsep.join([*prepend, *current_parts])


_augment_environment_path()


def require_tool(name: str) -> str:
    """Resolve a required CLI binary on PATH or fail the test (never skip)."""
    found = shutil.which(name)
    assert found is not None, (
        f"Required bioinformatics/workflow tool '{name}' not found on PATH: "
        f"{os.environ.get('PATH', '')}"
    )
    return found


@pytest.fixture()
def minimal_workspace(tmp_path: Path) -> dict[str, Any]:
    """Create a self-contained workspace with dummy BAMs, VCFs, BEDs, and FASTA."""
    bam_dir = tmp_path / "bams"
    bam_dir.mkdir()
    vcf_dir = tmp_path / "vcfs"
    (vcf_dir / "mutect2" / "filtered").mkdir(parents=True)
    (vcf_dir / "mutect2" / "merged").mkdir(parents=True)
    ref_dir = tmp_path / "ref"
    ref_dir.mkdir()
    out_dir = tmp_path / "results"

    fasta = ref_dir / "GRCh38.fa"
    fasta.write_text(">chr1\n" + ("ACGT" * 250) + "\n", encoding="utf-8")
    (ref_dir / "GRCh38.fa.fai").write_text("chr1\t1000\t6\t80\t81\n", encoding="utf-8")

    targets_bed = ref_dir / "targets.bed"
    targets_bed.write_text(
        "chr1\t100\t300\tGENE1_ex1\nchr1\t500\t700\tGENE2_ex1\n",
        encoding="utf-8",
    )

    access_bed = ref_dir / "access.bed"
    access_bed.write_text("chr1\t0\t1000\n", encoding="utf-8")

    blacklist_bed = ref_dir / "blacklist.bed"
    blacklist_bed.write_text("chr1\t250\t350\n", encoding="utf-8")

    ref_flat = ref_dir / "refFlat.txt"
    ref_flat.write_text(
        "GENE1\tNM_0001\tchr1\t+\t100\t300\t100\t300\t1\t100,\t300,\n",
        encoding="utf-8",
    )

    for bam_name in ("T1", "T2_R1", "T2_R2", "N1", "N2", "TO1", "G1"):
        (bam_dir / f"{bam_name}.merged.dedup.bqsr.bam").write_bytes(b"BAM\x01")
        (bam_dir / f"{bam_name}.merged.dedup.bqsr.bai").write_bytes(b"BAI\x01")

    for sname in ("PAT1", "PAT2", "PAT_TO"):
        vcf_file = vcf_dir / "mutect2" / "filtered" / f"{sname}.filtered.vcf.gz"
        vcf_file.write_bytes(b"\x1f\x8b")
        (vcf_dir / "mutect2" / "filtered" / f"{sname}.filtered.vcf.gz.tbi").write_bytes(b"TBI")
        (vcf_dir / "mutect2" / "merged" / f"{sname}.vcf.gz.stats").write_text(
            "stats\n", encoding="utf-8"
        )

    samples_tsv = tmp_path / "samples.tsv"
    samples_df = pd.DataFrame(
        [
            {
                "sample": "PAT1",
                "tumor_bam": "T1",
                "normal_bam": "N1",
                "analysis_type": "tumor_normal",
                "in_pon": "yes",
                "sex": "female",
            },
            {
                "sample": "PAT2",
                "tumor_bam": "T2_R1;T2_R2",
                "normal_bam": "N2",
                "analysis_type": "tumor_normal",
                "in_pon": "yes",
                "sex": "male",
            },
            {
                "sample": "PAT_TO",
                "tumor_bam": "TO1",
                "normal_bam": ".",
                "analysis_type": "tumor_only",
                "in_pon": "no",
                "sex": "female",
            },
            {
                "sample": "NORM_ONLY",
                "tumor_bam": "G1",
                "normal_bam": ".",
                "analysis_type": "germline",
                "in_pon": "yes",
                "sex": "male",
            },
        ]
    )
    samples_df.to_csv(samples_tsv, sep="\t", index=False)

    config_dict: dict[str, Any] = {
        "ref": {
            "genome": str(fasta),
            "build": "GRCh38",
            "targets_bed": str(targets_bed),
            "baits_bed": "",
            "access_bed": str(access_bed),
            "blacklist_bed": str(blacklist_bed),
            "ref_flat": str(ref_flat),
            "mappability_bw": "",
            "reptiming_bw": "",
        },
        "paths": {
            "samples": str(samples_tsv),
            "bam_folder": str(bam_dir),
            "vcf_folder": str(vcf_dir),
            "output_folder": str(out_dir),
            "log_subdir": "logs",
        },
        "bam": {
            "file_extension": ".merged.dedup.bqsr.bam",
        },
        "cnvkit": {
            "reference_cnn": "",
            "autobin": False,
            "target_avg_size": 200,
            "antitarget_avg_size": 150000,
            "antitarget_min_size": 500,
            "access_min_gap_size": 5000,
            "pon_bivar_threshold": 0.3,
            "male_reference": False,
            "segment_method": "cbs",
            "segment_threshold": 0.0001,
            "drop_low_coverage": True,
            "call_method": "clonal",
            "call_filter": "ci",
            "min_variant_depth": 20,
            "genes_of_interest": ["TP53", "EGFR"],
        },
        "purecn": {
            "enabled": True,
            "mode": "cnvkit",
            "genome": "hg38",
            "assay": "exome",
            "fun_segmentation": "PSCBS",
            "normal_panel_vcf": "",
            "normaldb": "",
            "mapping_bias": "",
            "snp_blacklist": str(blacklist_bed),
            "seed": 123,
            "postoptimize": True,
            "extra": "",
        },
        "params": {
            "cnvkit": {
                "segment_method": "cbs",
                "segment_threshold": 0.0001,
                "access_extra": "",
                "target_extra": "",
                "antitarget_extra": "",
                "autobin_extra": "",
                "coverage_extra": "",
                "reference_extra": "",
                "fix_extra": "",
                "segment_extra": "",
                "segmetrics_extra": "",
                "call_extra": "",
                "bintest_extra": "",
                "scatter_extra": "",
                "diagram_extra": "",
                "heatmap_extra": "",
            },
            "purecn": {
                "intervalfile_extra": "",
                "coverage_extra": "",
                "normaldb_extra": "",
                "run_extra": "",
            },
        },
    }

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config_dict, sort_keys=False), encoding="utf-8")

    return {
        "tmp_path": tmp_path,
        "bam_dir": bam_dir,
        "vcf_dir": vcf_dir,
        "ref_dir": ref_dir,
        "out_dir": out_dir,
        "fasta": fasta,
        "targets_bed": targets_bed,
        "access_bed": access_bed,
        "blacklist_bed": blacklist_bed,
        "ref_flat": ref_flat,
        "samples_tsv": samples_tsv,
        "samples_df": samples_df,
        "config_dict": config_dict,
        "config_path": config_path,
    }
