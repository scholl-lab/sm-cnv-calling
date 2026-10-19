"""Schema validation tests for config.schema.yaml and samples.schema.yaml."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml
from snakemake.exceptions import WorkflowError as _SchemaError
from snakemake.utils import validate as snakemake_validate

REPO_ROOT = Path(__file__).resolve().parent.parent
SCHEMA_DIR = REPO_ROOT / "workflow" / "schemas"
CONFIG_SCHEMA = str(SCHEMA_DIR / "config.schema.yaml")
SAMPLES_SCHEMA = str(SCHEMA_DIR / "samples.schema.yaml")


class TestConfigSchema:
    """Tests for workflow/schemas/config.schema.yaml."""

    def test_repo_default_config_valid(self) -> None:
        """The shipped config/config.yaml passes schema validation with GRCh38/hg38 defaults."""
        cfg = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text(encoding="utf-8"))
        snakemake_validate(cfg, CONFIG_SCHEMA)
        assert cfg["ref"]["build"] == "GRCh38"
        assert cfg["purecn"]["genome"] == "hg38"

    def test_minimal_workspace_config_valid(self, minimal_workspace: dict[str, Any]) -> None:
        """Minimal workspace fixture config validates cleanly."""
        snakemake_validate(minimal_workspace["config_dict"], CONFIG_SCHEMA)

    def test_missing_required_top_level_keys_rejected(
        self, minimal_workspace: dict[str, Any]
    ) -> None:
        """Removing required 'paths', 'ref', or 'bam' top-level sections fails schema validation."""
        for required_key in ("paths", "ref", "bam"):
            cfg = copy.deepcopy(minimal_workspace["config_dict"])
            del cfg[required_key]
            with pytest.raises(_SchemaError):
                snakemake_validate(cfg, CONFIG_SCHEMA)

    def test_missing_required_ref_fields_rejected(self, minimal_workspace: dict[str, Any]) -> None:
        """Removing required ref.genome or ref.targets_bed fails schema validation."""
        for field in ("genome", "targets_bed"):
            cfg = copy.deepcopy(minimal_workspace["config_dict"])
            del cfg["ref"][field]
            with pytest.raises(_SchemaError):
                snakemake_validate(cfg, CONFIG_SCHEMA)

    @pytest.mark.parametrize("bad_build", ["hg18", "GRCh36", "mm10", "invalid"])
    def test_invalid_ref_build_enum_rejected(
        self, minimal_workspace: dict[str, Any], bad_build: str
    ) -> None:
        """Unsupported ref.build values are rejected by schema enum."""
        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["ref"]["build"] = bad_build
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, CONFIG_SCHEMA)

    @pytest.mark.parametrize("bad_gv", ["GRCh38_invalid", "hg20", "dm6"])
    def test_invalid_purecn_genome_rejected(
        self, minimal_workspace: dict[str, Any], bad_gv: str
    ) -> None:
        """Unsupported purecn.genome values are rejected by schema enum."""
        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["purecn"]["genome"] = bad_gv
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, CONFIG_SCHEMA)

    def test_invalid_cnvkit_segment_method_rejected(
        self, minimal_workspace: dict[str, Any]
    ) -> None:
        """Invalid cnvkit.segment_method enum is rejected."""
        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["cnvkit"]["segment_method"] = "exome_invalid"
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, CONFIG_SCHEMA)

    def test_invalid_purecn_mode_rejected(self, minimal_workspace: dict[str, Any]) -> None:
        """Invalid purecn.mode enum is rejected."""
        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["purecn"]["mode"] = "invalid_mode"
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, CONFIG_SCHEMA)

    def test_invalid_pon_bivar_threshold_rejected(self, minimal_workspace: dict[str, Any]) -> None:
        """Negative pon_bivar_threshold is rejected."""
        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["cnvkit"]["pon_bivar_threshold"] = -0.1
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, CONFIG_SCHEMA)

    def test_unknown_top_level_property_rejected(self, minimal_workspace: dict[str, Any]) -> None:
        """Unknown top-level keys are rejected (additionalProperties: false)."""
        cfg = copy.deepcopy(minimal_workspace["config_dict"])
        cfg["nonexistent_section"] = {"foo": "bar"}
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, CONFIG_SCHEMA)


class TestSamplesSchema:
    """Tests for workflow/schemas/samples.schema.yaml."""

    def test_repo_default_samples_valid(self) -> None:
        """The shipped config/samples.tsv passes schema validation."""
        df = pd.read_csv(
            REPO_ROOT / "config" / "samples.tsv",
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        snakemake_validate(df, SAMPLES_SCHEMA)

    def test_valid_modern_and_legacy_rows(self) -> None:
        """Modern (tumor_normal, tumor_only, germline) and legacy (TvsN, To) rows pass."""
        df = pd.DataFrame(
            [
                {
                    "sample": "S1",
                    "tumor_bam": "T1;T2",
                    "normal_bam": "N1",
                    "analysis_type": "tumor_normal",
                    "in_pon": "yes",
                    "sex": "female",
                },
                {
                    "sample": "S2",
                    "tumor_bam": "T3",
                    "normal_bam": ".",
                    "analysis_type": "tumor_only",
                    "in_pon": "no",
                    "sex": "M",
                },
                {
                    "sample": "S3",
                    "tumor_bam": "N3",
                    "normal_bam": ".",
                    "analysis_type": "germline",
                    "in_pon": "true",
                    "sex": "F",
                },
                {
                    "sample": "S4",
                    "tumor_bam": "T4",
                    "normal_bam": "N4",
                    "analysis_type": "TvsN",
                },
                {
                    "sample": "S5",
                    "tumor_bam": "T5",
                    "normal_bam": ".",
                    "analysis_type": "To",
                },
            ]
        )
        snakemake_validate(df, SAMPLES_SCHEMA)

    def test_missing_required_sample_column_rejected(self) -> None:
        """Missing required 'sample' or 'analysis_type' column fails schema validation."""
        df = pd.DataFrame([{"tumor_bam": "T1", "normal_bam": "N1"}])
        with pytest.raises(_SchemaError):
            snakemake_validate(df, SAMPLES_SCHEMA)

    def test_invalid_analysis_type_rejected(self) -> None:
        """Invalid analysis_type enum fails schema validation."""
        df = pd.DataFrame(
            [{"sample": "S1", "tumor_bam": "T1", "normal_bam": "N1", "analysis_type": "RNAseq"}]
        )
        with pytest.raises(_SchemaError):
            snakemake_validate(df, SAMPLES_SCHEMA)
