"""Snakemake dry-run (-n -p) and launcher tests for sm-cnv-calling."""

from __future__ import annotations

import copy
import os
import re
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml

from tests.conftest import require_tool

pytestmark = [pytest.mark.dryrun]

REPO_ROOT = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"


def _run_snakemake(
    config_path: Path,
    *extra_args: str,
    dry_run: bool = True,
    cwd: Path = REPO_ROOT,
) -> subprocess.CompletedProcess[str]:
    """Execute snakemake as a subprocess and return CompletedProcess."""
    snk = require_tool("snakemake")
    cmd = [
        snk,
        "-s",
        str(SNAKEFILE),
        "--configfile",
        str(config_path),
        "--directory",
        str(cwd),
    ]
    if dry_run:
        cmd.append("-n")
    cmd.extend(extra_args)
    return subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=90,
    )


class TestWorkflowDryRun:
    """Dry-run tests verifying DAG construction and command flags across modes."""

    def test_list_rules(self, minimal_workspace: dict[str, Any]) -> None:
        """All modular CNVkit, PureCN, and QC rules are registered in the workflow."""
        res = _run_snakemake(minimal_workspace["config_path"], "--list-rules", dry_run=False)
        assert res.returncode == 0, f"--list-rules failed:\n{res.stdout}\n{res.stderr}"
        expected_rules = [
            "cnvkit_access",
            "cnvkit_filter_targets",
            "cnvkit_target",
            "cnvkit_antitarget",
            "cnvkit_coverage_target",
            "cnvkit_coverage_antitarget",
            "cnvkit_sex",
            "resolve_unit_sex",
            "cnvkit_metrics_normals",
            "identify_clean_normals",
            "cnvkit_reference_pooled",
            "cnvkit_fix",
            "cnvkit_segment",
            "cnvkit_segmetrics",
            "cnvkit_bintest",
            "cnvkit_call",
            "cnvkit_export_seg",
            "cnvkit_export_bed",
            "cnvkit_export_vcf",
            "cnvkit_scatter_genome",
            "cnvkit_scatter_gene",
            "cnvkit_diagram",
            "cnvkit_heatmap_cohort",
            "purecn_intervals",
            "purecn_coverage",
            "purecn_normaldb",
            "purecn_run",
            "consolidate_purity",
            "create_summary_report",
        ]
        for rname in expected_rules:
            assert rname in res.stdout, f"Missing expected rule '{rname}' in --list-rules"

    def test_full_cohort_dryrun_cnvkit_and_purecn_commands(
        self, minimal_workspace: dict[str, Any]
    ) -> None:
        """Dry-run (-n -p) verifies target/blacklist, modular CNVkit, BAF CBS, PureCN, and exports."""
        res = _run_snakemake(minimal_workspace["config_path"], "-p")
        assert res.returncode == 0, (
            f"Dry-run failed:\nstdout:\n{res.stdout}\nstderr:\n{res.stderr}"
        )
        out = res.stdout

        # #1: Blacklist subtraction on targets and access BEDs
        assert "rule cnvkit_filter_targets:" in out
        assert "rule cnvkit_access:" in out
        assert "bedtools subtract" in out

        # #2: cnvkit target always uses --split --short-names and optional --annotate / --avg-size
        assert "rule cnvkit_target:" in out
        assert "--split" in out and "--short-names" in out
        assert "--annotate " in out
        assert "--avg-size 200" in out

        # #3: Per-BAM coverage, fix, allele-aware CBS segment, and bintest
        assert "rule cnvkit_coverage_target:" in out
        assert "rule cnvkit_coverage_antitarget:" in out
        assert "rule cnvkit_fix:" in out
        assert "rule cnvkit_segment:" in out
        assert "-m cbs" in out and "-t 0.0001" in out
        assert "--drop-low-coverage" in out
        assert "-v " in out and "-i T1 -n N1" in out
        assert "rule cnvkit_bintest:" in out

        # #5: Multi-region units PAT2__T2_R1 and PAT2__T2_R2 scheduled with -i T2_R1 -n N2
        assert "PAT2__T2_R1.cnr" in out
        assert "PAT2__T2_R2.cnr" in out
        assert "-i T2_R1 -n N2" in out
        assert "-i T2_R2 -n N2" in out

        # #6: Normal QC metrics -> identify_clean_normals -> cnvkit_reference_pooled
        assert "rule cnvkit_metrics_normals:" in out
        assert "rule identify_clean_normals:" in out
        assert "rule cnvkit_reference_pooled:" in out

        # #7: PureCN cnvkit mode runs for both tumor_normal and tumor_only units
        assert "rule purecn_run:" in out
        assert "--tumor " in out and "PAT1.cnr" in out
        assert "--seg-file " in out and "PAT1.seg" in out
        assert "PAT_TO.cnr" in out
        assert "--statsfile " in out
        assert "--snpblacklist " in out
        assert "--genome hg38" in out
        assert "--fun-segmentation PSCBS" in out

        # #8: segmetrics, call, exports (seg, bed, vcf), gene scatter, summary
        assert "rule cnvkit_segmetrics:" in out
        assert "--ci --pi --sem --iqr" in out
        assert "rule cnvkit_call:" in out
        assert "rule cnvkit_export_seg:" in out
        assert "rule cnvkit_export_bed:" in out
        assert "rule cnvkit_export_vcf:" in out
        assert "rule cnvkit_scatter_gene:" in out
        assert "rule create_summary_report:" in out

    def test_auto_access_bed_and_autobin_mode(self, minimal_workspace: dict[str, Any]) -> None:
        """When ref.access_bed is empty and cnvkit.autobin is true, cnvkit_access and cnvkit_autobin run."""
        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["ref"]["access_bed"] = ""
        cfg["cnvkit"]["autobin"] = True
        cfg_path = minimal_workspace["tmp_path"] / "config_autobin.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

        res = _run_snakemake(cfg_path, "-p")
        assert res.returncode == 0, f"Autobin dry-run failed:\n{res.stdout}\n{res.stderr}"
        assert "rule cnvkit_access:" in res.stdout
        assert "cnvkit.py access" in res.stdout
        assert "rule cnvkit_autobin:" in res.stdout
        assert "cnvkit.py autobin" in res.stdout
        assert "rule cnvkit_target:" not in res.stdout

    def test_zero_normal_cohort_uses_flat_reference(
        self, minimal_workspace: dict[str, Any]
    ) -> None:
        """Cohort with 0 PoN normal BAMs schedules cnvkit_reference_pooled in flat-reference mode."""
        samples_to = minimal_workspace["tmp_path"] / "samples_to.tsv"
        pd.DataFrame(
            [
                {
                    "sample": "PAT_TO",
                    "tumor_bam": "TO1",
                    "normal_bam": ".",
                    "analysis_type": "tumor_only",
                    "in_pon": "no",
                }
            ]
        ).to_csv(samples_to, sep="\t", index=False)

        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["paths"]["samples"] = str(samples_to)
        cfg_path = minimal_workspace["tmp_path"] / "config_to.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

        res = _run_snakemake(cfg_path, "-p")
        assert res.returncode == 0, f"0-normal dry-run failed:\n{res.stdout}\n{res.stderr}"
        assert "rule cnvkit_reference_pooled:" in res.stdout
        assert "rule cnvkit_metrics_normals:" not in res.stdout
        assert "rule identify_clean_normals:" not in res.stdout

    def test_prebuilt_reference_cnn_bypasses_pon_rules(
        self, minimal_workspace: dict[str, Any]
    ) -> None:
        """Setting cnvkit.reference_cnn bypasses PoN coverage and reference construction."""
        prebuilt = minimal_workspace["ref_dir"] / "prebuilt_ref.cnn"
        prebuilt.write_text("cnn\n", encoding="utf-8")

        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["cnvkit"]["reference_cnn"] = str(prebuilt)
        cfg_path = minimal_workspace["tmp_path"] / "config_prebuilt.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

        res = _run_snakemake(cfg_path, "-p")
        assert res.returncode == 0, (
            f"Prebuilt reference dry-run failed:\n{res.stdout}\n{res.stderr}"
        )
        assert "rule cnvkit_reference_pooled:" not in res.stdout
        assert "rule cnvkit_metrics_normals:" not in res.stdout
        assert str(prebuilt) in res.stdout

    def test_purecn_native_mode_schedules_intervalfile_coverage_normaldb(
        self, minimal_workspace: dict[str, Any]
    ) -> None:
        """purecn.mode='native' schedules purecn_intervals, purecn_coverage, and purecn_normaldb."""
        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["purecn"]["mode"] = "native"
        cfg_path = minimal_workspace["tmp_path"] / "config_purecn_native.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

        res = _run_snakemake(cfg_path, "-p")
        assert res.returncode == 0, f"PureCN native dry-run failed:\n{res.stdout}\n{res.stderr}"
        out = res.stdout
        assert "rule purecn_intervals:" in out
        assert "IntervalFile.R" in out and "--offtarget" in out
        assert "rule purecn_coverage:" in out
        assert "rule purecn_normaldb:" in out
        assert "NormalDB.R" in out and "--coveragefiles" in out
        assert "--normaldb " in out and "--intervals " in out

    def test_default_profile_and_bih_profile_layering(
        self, minimal_workspace: dict[str, Any]
    ) -> None:
        """Layering profiles/default and profiles/bih preserves slurm_account and keeps runtimes <= 480m."""
        bih_dir = minimal_workspace["tmp_path"] / "bih_profile"
        bih_dir.mkdir()
        bih_lines = [
            line
            for line in (REPO_ROOT / "profiles" / "bih" / "config.yaml")
            .read_text(encoding="utf-8")
            .splitlines()
            if not line.strip().startswith("executor:")
        ]
        (bih_dir / "config.yaml").write_text("\n".join(bih_lines) + "\n", encoding="utf-8")

        res = _run_snakemake(
            minimal_workspace["config_path"],
            "--workflow-profile",
            str(REPO_ROOT / "profiles" / "default"),
            "--profile",
            str(bih_dir),
        )
        assert res.returncode == 0, f"Profile dry-run failed:\n{res.stdout}\n{res.stderr}"
        res_lines = [
            line.strip()
            for line in res.stdout.splitlines()
            if line.strip().startswith("resources:")
        ]
        assert res_lines, "Expected resources: lines in dry-run output"
        for line in res_lines:
            assert "slurm_account=hpc-ag-scholl" in line, f"Missing slurm_account in {line}"
            assert "slurm_partition=standard" in line, f"Missing slurm_partition in {line}"
            m = re.search(r"\bruntime=(\d+)\b", line)
            assert m is not None, f"Missing numeric runtime in {line}"
            assert int(m.group(1)) <= 480, f"runtime {m.group(1)} > 480 in {line}"

    def test_run_snakemake_launcher_script(self, minimal_workspace: dict[str, Any]) -> None:
        """scripts/run_snakemake.sh resolves REPO_ROOT, profiles, and forwards arguments cleanly."""
        stage_dir = minimal_workspace["tmp_path"] / "stage"
        stage_dir.mkdir()
        bin_dir = minimal_workspace["tmp_path"] / "bin"
        bin_dir.mkdir()
        args_file = minimal_workspace["tmp_path"] / "snakemake_args.txt"
        (bin_dir / "snakemake").write_text(
            f'#!/bin/bash\nprintf "%s\\n" "$@" > "{args_file}"\n',
            encoding="utf-8",
        )
        (bin_dir / "snakemake").chmod(0o755)

        launcher = REPO_ROOT / "scripts" / "run_snakemake.sh"
        res = subprocess.run(
            [
                "bash",
                str(launcher),
                str(SNAKEFILE),
                str(minimal_workspace["config_path"]),
                "--dry-run",
            ],
            cwd=stage_dir,
            env={**os.environ, "PATH": f"{bin_dir}:/usr/bin:/bin"},
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert res.returncode == 0, f"Launcher failed:\n{res.stdout}\n{res.stderr}"
        logged = args_file.read_text(encoding="utf-8").splitlines()
        assert str(SNAKEFILE) in logged
        assert str(minimal_workspace["config_path"]) in logged
        assert str(REPO_ROOT / "profiles" / "default") in logged
        assert "--dry-run" in logged
