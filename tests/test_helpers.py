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


class TestBedtoolsBlacklistSubtraction:
    """Tests for Issue #1: real bedtools subtract execution on target and access BEDs."""

    def test_bedtools_subtract_excises_blacklist_intervals(self, tmp_path: Path) -> None:
        """Run real bedtools subtract on targets.bed and blacklist.bed."""
        bedtools = require_tool("bedtools")
        targets = tmp_path / "targets.bed"
        blacklist = tmp_path / "blacklist.bed"
        clean = tmp_path / "targets.clean.bed"

        targets.write_text(
            "chr1\t100\t300\tEXON1\nchr1\t500\t700\tEXON2\n",
            encoding="utf-8",
        )
        blacklist.write_text("chr1\t250\t350\n", encoding="utf-8")

        res = subprocess.run(
            [bedtools, "subtract", "-a", str(targets), "-b", str(blacklist)],
            check=True,
            capture_output=True,
            text=True,
        )
        clean.write_text(res.stdout, encoding="utf-8")
        lines = [line.strip() for line in clean.read_text(encoding="utf-8").splitlines() if line]
        assert lines == [
            "chr1\t100\t250\tEXON1",
            "chr1\t500\t700\tEXON2",
        ]

    def test_samtools_bgzip_tabix_bcftools_execution(self, tmp_path: Path) -> None:
        """Run real samtools faidx, bgzip, tabix -p vcf, and bcftools view on synthetic CNV VCF."""
        samtools = require_tool("samtools")
        bgzip = require_tool("bgzip")
        tabix = require_tool("tabix")
        bcftools = require_tool("bcftools")

        fasta = tmp_path / "ref.fa"
        fasta.write_text(">chr1\n" + ("ACGT" * 50) + "\n", encoding="utf-8")
        subprocess.run([samtools, "faidx", str(fasta)], check=True, capture_output=True)
        assert (tmp_path / "ref.fa.fai").is_file()

        raw_vcf = tmp_path / "sample.cnv.vcf"
        raw_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=200>\n"
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPAT1\n"
            "chr1\t10\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=150\tGT\t0/1\n",
            encoding="utf-8",
        )
        vcf_gz = tmp_path / "sample.cnv.vcf.gz"
        with vcf_gz.open("wb") as out_fh:
            subprocess.run([bgzip, "-c", str(raw_vcf)], check=True, stdout=out_fh)
        subprocess.run([tabix, "-p", "vcf", str(vcf_gz)], check=True, capture_output=True)
        assert (tmp_path / "sample.cnv.vcf.gz.tbi").is_file()

        view_res = subprocess.run(
            [bcftools, "view", "-H", str(vcf_gz), "chr1:1-200"],
            check=True,
            capture_output=True,
            text=True,
        )
        assert "SVTYPE=DUP;END=150" in view_res.stdout


class TestNormalPoNFiltering:
    """Tests for Issue #6: normal QC metrics filtering and path resolution."""

    def test_identify_clean_normals_resolves_target_and_antitarget_cnn(
        self, tmp_path: Path
    ) -> None:
        """identify_clean_normals strips suffixes and outputs verified target+antitarget .cnn paths."""
        cov_dir = tmp_path / "coverage"
        cov_dir.mkdir()
        for bn in ("N1", "N2", "N3"):
            (cov_dir / f"{bn}.targetcoverage.cnn").write_text("cov\n", encoding="utf-8")
            (cov_dir / f"{bn}.antitargetcoverage.cnn").write_text("cov\n", encoding="utf-8")

        metrics_tsv = tmp_path / "normals_metrics.tsv"
        pd.DataFrame(
            [
                {"sample": str(cov_dir / "N1.targetcoverage.cnn"), "bivar": 0.12, "mad": 0.1},
                {"sample": "N2.cnn", "bivar": 0.45, "mad": 0.4},  # noisy -> excluded
                {"sample": "N3", "bivar": 0.25, "mad": 0.2},
            ]
        ).to_csv(metrics_tsv, sep="\t", index=False)

        out_list = tmp_path / "clean_normals.txt"
        paths = identify_clean_normals(metrics_tsv, cov_dir, out_list, threshold=0.3)

        assert paths == [
            str(cov_dir / "N1.targetcoverage.cnn"),
            str(cov_dir / "N1.antitargetcoverage.cnn"),
            str(cov_dir / "N3.targetcoverage.cnn"),
            str(cov_dir / "N3.antitargetcoverage.cnn"),
        ]
        assert out_list.read_text(encoding="utf-8").splitlines() == paths

    def test_identify_clean_normals_empty_or_all_noisy(self, tmp_path: Path) -> None:
        """Empty metrics file or all-noisy normals produce an empty clean_normals.txt without crashing."""
        cov_dir = tmp_path / "coverage"
        cov_dir.mkdir()
        empty_metrics = tmp_path / "empty.tsv"
        empty_metrics.write_text("", encoding="utf-8")
        out_list = tmp_path / "clean_normals.txt"

        assert identify_clean_normals(empty_metrics, cov_dir, out_list, threshold=0.3) == []
        assert out_list.read_text(encoding="utf-8") == ""

        noisy_metrics = tmp_path / "noisy.tsv"
        pd.DataFrame([{"sample": "N1", "bivar": 0.85}]).to_csv(
            noisy_metrics, sep="\t", index=False
        )
        assert identify_clean_normals(noisy_metrics, cov_dir, out_list, threshold=0.3) == []
        assert out_list.read_text(encoding="utf-8") == ""


class TestPureCNConsolidationAndCallArgs:
    """Tests for Issue #7 & #8: purity/ploidy consolidation, sex inference, and cnvkit call flags."""

    def test_consolidate_purity_ploidy_from_purecn_csv(self, tmp_path: Path) -> None:
        """Valid PureCN CSV extracts top-ranked Purity and Ploidy."""
        csv_file = tmp_path / "PAT1.csv"
        pd.DataFrame(
            [{"Sampleid": "PAT1", "Purity": 0.64, "Ploidy": 3.12, "Flagged": False}]
        ).to_csv(csv_file, index=False)
        status_file = tmp_path / "PAT1.status"
        status_file.write_text("succeeded\n", encoding="utf-8")
        out_tsv = tmp_path / "PAT1.purity_ploidy.tsv"

        rec = consolidate_purity_ploidy(
            unit="PAT1",
            sample="PAT1",
            purecn_csv=csv_file,
            purecn_status_file=status_file,
            sample_purity="",
            sample_ploidy="",
            output_file=out_tsv,
            purecn_attempted=True,
        )
        assert rec["purity"] == 0.64
        assert rec["ploidy"] == 3.12
        assert rec["source"] == "PureCN"
        assert rec["status"] == "succeeded"

        call_args = build_cnvkit_call_args(
            out_tsv,
            sex="female",
            male_reference=True,
            call_method="clonal",
            call_filter="ci",
            vcf_args="-v pat1.vcf.gz -i T1 -n N1",
        )
        assert "--filter ci" in call_args
        assert "-m clonal --purity 0.64 --ploidy 3" in call_args
        assert "-y" in call_args
        assert "-x female" in call_args
        assert "-v pat1.vcf.gz -i T1 -n N1" in call_args

    def test_consolidate_purity_fallback_to_threshold_when_unestimated(
        self, tmp_path: Path
    ) -> None:
        """When PureCN fails and no sample_sheet purity is given, cnvkit call falls back to -m threshold."""
        csv_file = tmp_path / "PAT_FAIL.csv"
        csv_file.write_text("", encoding="utf-8")
        status_file = tmp_path / "PAT_FAIL.status"
        status_file.write_text("failed\n", encoding="utf-8")
        out_tsv = tmp_path / "PAT_FAIL.purity_ploidy.tsv"

        rec = consolidate_purity_ploidy(
            unit="PAT_FAIL",
            sample="PAT_FAIL",
            purecn_csv=csv_file,
            purecn_status_file=status_file,
            sample_purity="",
            sample_ploidy="",
            output_file=out_tsv,
            purecn_attempted=True,
        )
        assert rec["purity"] == ""
        assert rec["ploidy"] == ""
        assert rec["source"] == "unestimated"
        assert rec["status"] == "failed"

        call_args = build_cnvkit_call_args(
            out_tsv,
            sex="male",
            male_reference=False,
            call_method="clonal",
            call_filter="ci",
        )
        assert "-m threshold" in call_args
        assert "--purity" not in call_args
        assert "-x male" in call_args

    def test_sex_parsing_and_resolution(self, tmp_path: Path) -> None:
        """parse_cnvkit_sex handles sample sheet overrides and cnvkit sex TSV output."""
        assert normalize_sex("F") == "female"
        assert normalize_sex("M") == "male"
        with pytest.raises(ValueError, match="Invalid sex"):
            normalize_sex("unknown_sex")

        sex_out = "sample\tsex\tchrX_ratio\tchrY_ratio\nT1.cnr\tMale\t-0.95\t-0.10\n"
        inferred_file = tmp_path / "T1.sex.tsv"
        inferred_file.write_text(sex_out, encoding="utf-8")
        assert parse_cnvkit_sex(inferred_file, "") == "male"
        assert parse_cnvkit_sex(inferred_file, "female") == "female"


class TestSummaryReportAndConfigGenerator:
    """Tests for pipeline summary generation and scripts/generate_config.py."""
