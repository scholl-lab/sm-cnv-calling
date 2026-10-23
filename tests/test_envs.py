"""Static checks for the pinned CNVkit segmentation environment."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SHA = "7539ce336725351706f9adfe32b889b05252b3e7"


def test_cnvkit_segment_rule_uses_pinned_env() -> None:
    smk = (ROOT / "workflow/rules/cnvkit.smk").read_text()
    block = re.search(r"rule cnvkit_segment:.*?shell:", smk, re.S)
    assert block is not None
    assert '"../envs/cnvkit_segment.yaml"' in block.group(0)


def test_cnvkit_segment_post_deploy_pins_sha_without_deps() -> None:
    envs = ROOT / "workflow/envs"
    script = envs / "cnvkit_segment.post-deploy.sh"
    text = script.read_text()
    assert f'pip install --no-deps "git+https://github.com/etal/cnvkit@{SHA}"' in text
    assert "0.9.15.dev0" in text
    assert (envs / "cnvkit_segment.yaml").read_text() == (envs / "cnvkit.yaml").read_text()
