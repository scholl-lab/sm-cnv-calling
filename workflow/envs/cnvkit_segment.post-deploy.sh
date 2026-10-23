#!/usr/bin/env bash
# CNVkit 0.9.14 `segment -m cbs -v VCF` concatenates per-segment BAF results without
# ignore_index, which corrupts segment coordinates; fixed unreleased in etal/cnvkit#1125.
# ponytail: drop this env and use cnvkit.yaml again once a CNVkit release > 0.9.14 contains #1125.
set -euo pipefail
pip install --no-deps "git+https://github.com/etal/cnvkit@7539ce336725351706f9adfe32b889b05252b3e7"
cnvkit.py version | grep -q "0.9.15.dev0" || { echo "cnvkit_segment env: expected CNVkit 0.9.15.dev0" >&2; exit 1; }
