"""Unit and tool-execution tests for workflow/rules/helpers.py and scripts/generate_config.py."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml

from scripts.generate_config import (
    build_sample_rows,
    discover_bam_files,
    render_config_yaml,
)
from tests.conftest import require_tool
from workflow.rules.helpers import (
    build_cnvkit_call_args,
    build_cnvkit_vcf_args,
    build_pipeline_summary,
    consolidate_purity_ploidy,
    get_cnv_units,
    get_normal_bam_basenames,
    get_tumor_samples,
    get_unique_bam_basenames,
    identify_clean_normals,
    normalize_analysis_type,
    normalize_sex,
    parse_cnvkit_sex,
    resolve_bam_path,
    resolve_sample_vcf,
    resolve_sample_vcf_stats,
    strip_bam_extension,
    validate_samples_df,
)


class TestSampleSheetAndBamResolution:
    """Tests for Issue #5: sm-calling sample sheet alignment and 97-BAM cohort support."""

    def test_strip_and_resolve_bam_paths(self, tmp_path: Path) -> None:
        """Basename and explicit path BAM resolution matches sm-calling conventions."""
        ext = ".merged.dedup.bqsr.bam"
        assert strip_bam_extension("T1.merged.dedup.bqsr.bam", ext) == "T1"
        assert strip_bam_extension("/path/to/T1.bam", ext) == "T1"
        assert strip_bam_extension("T1", ext) == "T1"

        bam_dir = tmp_path / "bams"
        bam_dir.mkdir()
        resolved_bam = resolve_bam_path("T1", str(bam_dir), ext)
        assert resolved_bam == str(bam_dir / "T1.merged.dedup.bqsr.bam")

        explicit_bam = str(tmp_path / "custom" / "S9.bam")
        assert resolve_bam_path(explicit_bam, str(bam_dir), ext) == explicit_bam

    def test_normalize_analysis_type(self) -> None:
        """Modern and legacy analysis_type values normalize accurately."""
        assert normalize_analysis_type("tumor_normal") == "tumor_normal"
        assert normalize_analysis_type("TvsN") == "tumor_normal"
        assert normalize_analysis_type("tumor_only") == "tumor_only"
        assert normalize_analysis_type("To") == "tumor_only"
        assert normalize_analysis_type("germline") == "germline"
        with pytest.raises(ValueError, match="Invalid analysis_type"):
            normalize_analysis_type("unknown_mode")

    def test_multi_region_and_97_bam_cohort_expansion(self) -> None:
        """Simulate a 97-BAM cohort with multi-region tumors, shared normals, tumor-only, and germline PoN."""
        rows: list[dict[str, Any]] = []
        # 15 multi-region tumor_normal patients: each has 2 tumor BAMs + 1 normal BAM = 45 BAMs
        for i in range(1, 16):
            rows.append(
                {
                    "sample": f"MR_PAT_{i:02d}",
                    "tumor_bam": f"MR_T_{i:02d}A;MR_T_{i:02d}B",
                    "normal_bam": f"MR_N_{i:02d}",
                    "analysis_type": "tumor_normal",
                    "in_pon": "yes",
                }
            )
        # 20 single-region tumor_normal patients: 20 tumor BAMs + 20 normal BAMs = 40 BAMs
        for i in range(1, 21):
            rows.append(
                {
                    "sample": f"SR_PAT_{i:02d}",
                    "tumor_bam": f"SR_T_{i:02d}",
                    "normal_bam": f"SR_N_{i:02d}",
                    "analysis_type": "tumor_normal",
                    "in_pon": "yes",
                }
            )
        # 7 tumor_only patients: 7 tumor BAMs
        for i in range(1, 8):
            rows.append(
                {
                    "sample": f"TO_PAT_{i:02d}",
                    "tumor_bam": f"TO_T_{i:02d}",
                    "normal_bam": ".",
                    "analysis_type": "tumor_only",
                    "in_pon": "no",
                }
            )
        # 5 germline PoN-only normals: 5 normal BAMs -> total unique BAMs = 45 + 40 + 7 + 5 = 97 BAMs
        for i in range(1, 6):
            rows.append(
                {
                    "sample": f"GL_NORM_{i:02d}",
                    "tumor_bam": f"GL_N_{i:02d}",
                    "normal_bam": ".",
                    "analysis_type": "germline",
                    "in_pon": "yes",
                }
            )

        df = pd.DataFrame(rows)
        validate_samples_df(df)

        all_bams = get_unique_bam_basenames(df)
        somatic_samples = get_tumor_samples(df)
        pon_bams = get_normal_bam_basenames(df)
        units = get_cnv_units(df)

        assert len(all_bams) == 97
        assert len(somatic_samples) == 15 + 20 + 7  # 42 somatic samples
        assert len(pon_bams) == 15 + 20 + 5  # 40 PoN normal BAMs
        assert len(units) == 15 * 2 + 20 + 7  # 57 somatic CNV units (germline excluded)

        # Multi-region unit naming and VCF sample ID defaults
        u_mr = units["MR_PAT_01__MR_T_01A"]
        assert u_mr["sample"] == "MR_PAT_01"
        assert u_mr["tumor_bam"] == "MR_T_01A"
        assert u_mr["normal_bam"] == "MR_N_01"
        assert u_mr["tumor_vcf_id"] == "MR_T_01A"
        assert u_mr["normal_vcf_id"] == "MR_N_01"

        # Single-region unit keeps sample name as unit ID
        u_sr = units["SR_PAT_01"]
        assert u_sr["tumor_bam"] == "SR_T_01"
        assert u_sr["normal_bam"] == "SR_N_01"

        # Tumor-only unit has empty normal_vcf_id
        u_to = units["TO_PAT_01"]
        assert u_to["analysis_type"] == "tumor_only"
        assert u_to["normal_bam"] == ""
        assert u_to["normal_vcf_id"] == ""

    def test_validate_samples_errors(self) -> None:
        """Duplicate sample IDs or missing BAMs raise clear ValueErrors."""
        dup_df = pd.DataFrame(
            [
                {
                    "sample": "S1",
                    "tumor_bam": "T1",
                    "normal_bam": "N1",
                    "analysis_type": "tumor_normal",
                },
                {
                    "sample": "S1",
                    "tumor_bam": "T2",
                    "normal_bam": "N2",
                    "analysis_type": "tumor_normal",
                },
            ]
        )
        with pytest.raises(ValueError, match="Duplicate sample IDs"):
            validate_samples_df(dup_df)

        missing_normal_df = pd.DataFrame(
            [
                {
                    "sample": "S1",
                    "tumor_bam": "T1",
                    "normal_bam": ".",
                    "analysis_type": "tumor_normal",
                }
            ]
        )
        with pytest.raises(ValueError, match="missing normal_bam"):
            validate_samples_df(missing_normal_df)

    def test_vcf_and_stats_resolution(self, minimal_workspace: dict[str, Any]) -> None:
        """VCF and Mutect2 .vcf.gz.stats resolve from paths.vcf_folder or explicit columns."""
        df = minimal_workspace["samples_df"]
        vcf_dir = str(minimal_workspace["vcf_dir"])
        assert resolve_sample_vcf(df, "PAT1", vcf_dir).endswith(
            "mutect2/filtered/PAT1.filtered.vcf.gz"
        )
        assert resolve_sample_vcf_stats(df, "PAT1", vcf_dir).endswith(
            "mutect2/merged/PAT1.vcf.gz.stats"
        )
        assert resolve_sample_vcf(df, "PAT1", "") == ""

        # Explicit VCF column overrides vcf_folder
        df_explicit = df.copy()
        df_explicit["vcf"] = ["/custom/PAT1.vcf.gz", "", "", ""]
        df_explicit["vcf_stats"] = ["/custom/PAT1.vcf.gz.stats", "", "", ""]
        assert resolve_sample_vcf(df_explicit, "PAT1", vcf_dir) == "/custom/PAT1.vcf.gz"
        assert (
            resolve_sample_vcf_stats(df_explicit, "PAT1", vcf_dir) == "/custom/PAT1.vcf.gz.stats"
        )

    def test_build_cnvkit_vcf_args(self) -> None:
        """build_cnvkit_vcf_args formats -v, -i, -n, -m flags properly."""
        assert build_cnvkit_vcf_args("") == ""
        assert (
            build_cnvkit_vcf_args("pat1.vcf.gz", "T1", "N1", 20)
            == "-v pat1.vcf.gz -i T1 -n N1 -m 20"
        )
        assert build_cnvkit_vcf_args("to1.vcf.gz", "TO1", "", 15) == "-v to1.vcf.gz -i TO1 -m 15"
