"""End-to-end integration test executing the real Snakemake workflow with CNVkit 0.9.14."""

from __future__ import annotations

import copy
import random
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

CONTIGS: dict[str, int] = {
    "chr1": 16_000,
    "chr2": 16_000,
    "chr3": 24_000,
    "chrX": 16_000,
}
GAIN_REGION = ("chr1", 1_000, 13_000)
LOSS_REGION = ("chr2", 1_000, 13_000)


def _region_multiplier(chrom: str, pos: int, *, is_tumor: bool) -> float:
    """Return copy-ratio depth multiplier for a genomic coordinate."""
    if not is_tumor:
        return 1.0
    if chrom == GAIN_REGION[0] and GAIN_REGION[1] <= pos < GAIN_REGION[2]:
        return 2.6
    if chrom == LOSS_REGION[0] and LOSS_REGION[1] <= pos < LOSS_REGION[2]:
        return 0.35
    return 1.0


def _write_synthetic_bam(
    samtools: str,
    bam_path: Path,
    seqs: dict[str, str],
    targets: list[tuple[str, int, int, str]],
    *,
    seed: int,
    is_tumor: bool,
) -> None:
    """Generate a coordinate-sorted, indexed BAM with deterministic on/off-target depth."""
    rng = random.Random(seed)
    read_len = 200
    qual = "I" * read_len

    lines: list[str] = ["@HD\tVN:1.6\tSO:coordinate"]
    for chrom, length in CONTIGS.items():
        lines.append(f"@SQ\tSN:{chrom}\tLN:{length}")

    target_by_pos: dict[tuple[str, int], int] = {}
    for chrom, start, _end, _gene in targets:
        mult = _region_multiplier(chrom, start, is_tumor=is_tumor)
        jitter = rng.choice([-2, -1, 0, 1, 2])
        target_by_pos[(chrom, start)] = max(2, round(40 * mult) + jitter)

    qname_idx = 0
    for chrom, length in CONTIGS.items():
        seq = seqs[chrom]
        for start in range(0, length - read_len + 1, read_len):
            mult = _region_multiplier(chrom, start, is_tumor=is_tumor)
            n_reads = target_by_pos.get((chrom, start), max(1, round(6 * mult)))
            read_seq = seq[start : start + read_len]
            pos_1based = start + 1
            for _ in range(n_reads):
                qname_idx += 1
                lines.append(
                    f"r{qname_idx}\t0\t{chrom}\t{pos_1based}\t60\t{read_len}M\t*\t0\t0\t{read_seq}\t{qual}"
                )

    sam_text = "\n".join(lines) + "\n"
    subprocess.run(
        [samtools, "view", "-bS", "-", "-o", str(bam_path)],
        input=sam_text.encode("ascii"),
        check=True,
        capture_output=True,
    )
    subprocess.run([samtools, "index", str(bam_path)], check=True, capture_output=True)


def _write_germline_baf_vcf(
    bgzip: str,
    tabix: str,
    vcf_gz_path: Path,
    seqs: dict[str, str],
    targets: list[tuple[str, int, int, str]],
) -> None:
    """Create a bgzipped, tabix-indexed paired VCF with heterozygous SNPs for BAF."""
    alt_map = {"A": "G", "G": "A", "C": "T", "T": "C"}
    header = [
        "##fileformat=VCFv4.2",
        *(f"##contig=<ID={c},length={ln}>" for c, ln in CONTIGS.items()),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT1\tN1",
    ]
    records: list[str] = []
    for chrom, start, _end, _gene in targets[::2]:
        pos0 = start + 100
        pos1 = pos0 + 1
        ref_base = seqs[chrom][pos0]
        alt_base = alt_map[ref_base]
        n_gt = "0/1:30,30:60"
        if chrom == GAIN_REGION[0] and GAIN_REGION[1] <= start < GAIN_REGION[2]:
            t_gt = "0/1:30,90:120"
        elif chrom == LOSS_REGION[0] and LOSS_REGION[1] <= start < LOSS_REGION[2]:
            t_gt = "0/1:36,4:40"
        else:
            t_gt = "0/1:30,30:60"
        records.append(
            f"{chrom}\t{pos1}\t.\t{ref_base}\t{alt_base}\t60\tPASS\t.\tGT:AD:DP\t{t_gt}\t{n_gt}"
        )

    vcf_text = "\n".join([*header, *records]) + "\n"
    vcf_gz_path.parent.mkdir(parents=True, exist_ok=True)
    with vcf_gz_path.open("wb") as fh:
        subprocess.run(
            [bgzip, "-c"],
            input=vcf_text.encode("ascii"),
            stdout=fh,
            check=True,
        )
    subprocess.run([tabix, "-p", "vcf", str(vcf_gz_path)], check=True, capture_output=True)


def test_end_to_end_cnvkit_workflow(tmp_path: Path) -> None:
    """Run the full Snakemake CNVkit workflow and flat-reference fallback on synthetic data."""
    snakemake = require_tool("snakemake")
    cnvkit = require_tool("cnvkit.py")
    samtools = require_tool("samtools")
    require_tool("bedtools")
    bcftools = require_tool("bcftools")
    bgzip = require_tool("bgzip")
    tabix = require_tool("tabix")
    require_tool("Rscript")

    ver_res = subprocess.run([cnvkit, "version"], check=True, capture_output=True, text=True)
    assert ver_res.stdout.strip() == "0.9.14", f"Expected cnvkit 0.9.14, got {ver_res.stdout!r}"

    ref_dir = tmp_path / "ref"
    bam_dir = tmp_path / "bams"
    vcf_dir = tmp_path / "vcfs"
    out_dir = tmp_path / "results"
    for d in (ref_dir, bam_dir, vcf_dir):
        d.mkdir(parents=True)

    # 1. Tiny deterministic reference FASTA (>=3 contigs including chrX)
    rng = random.Random(42)
    seqs = {chrom: "".join(rng.choices("ACGT", k=length)) for chrom, length in CONTIGS.items()}
    fasta = ref_dir / "GRCh38.fa"
    fasta.write_text(
        "".join(f">{chrom}\n{seq}\n" for chrom, seq in seqs.items()),
        encoding="utf-8",
    )
    subprocess.run([samtools, "faidx", str(fasta)], check=True, capture_output=True)

    # 2. Capture targets BED (all start positions aligned to 200 bp windows), blacklist, refFlat
    targets: list[tuple[str, int, int, str]] = []
    for i in range(12):
        s = 1_000 + i * 1_000
        targets.append(("chr1", s, s + 200, "EGFR"))
    for i in range(12):
        s = 1_000 + i * 1_000
        targets.append(("chr2", s, s + 200, "TP53"))
    for i in range(18):
        s = 1_000 + i * 1_000
        targets.append(("chr3", s, s + 200, "BRCA1"))
    for i in range(12):
        s = 1_000 + i * 1_000
        targets.append(("chrX", s, s + 200, "AR"))

    targets_bed = ref_dir / "targets.bed"
    targets_bed.write_text(
        "".join(f"{c}\t{s}\t{e}\t{g}\n" for c, s, e, g in targets),
        encoding="utf-8",
    )

    blacklist_bed = ref_dir / "blacklist.bed"
    blacklist_bed.write_text("chr3\t23600\t23800\n", encoding="utf-8")

    ref_flat = ref_dir / "refFlat.txt"
    ref_flat.write_text(
        "EGFR\tNM_005228\tchr1\t+\t1000\t12200\t1000\t12200\t1\t1000,\t12200,\n"
        "TP53\tNM_000546\tchr2\t+\t1000\t12200\t1000\t12200\t1\t1000,\t12200,\n"
        "BRCA1\tNM_007294\tchr3\t+\t1000\t18200\t1000\t18200\t1\t1000,\t18200,\n"
        "AR\tNM_000044\tchrX\t+\t1000\t12200\t1000\t12200\t1\t1000,\t12200,\n",
        encoding="utf-8",
    )

    # 3. Generate 3 normal BAMs (N1, N2, N3) and 1 tumor BAM (T1 with gain + loss)
    ext = ".merged.dedup.bqsr.bam"
    for idx, norm_name in enumerate(("N1", "N2", "N3"), start=1):
        _write_synthetic_bam(
            samtools,
            bam_dir / f"{norm_name}{ext}",
            seqs,
            targets,
            seed=100 + idx,
            is_tumor=False,
        )
    _write_synthetic_bam(
        samtools,
        bam_dir / f"T1{ext}",
        seqs,
        targets,
        seed=999,
        is_tumor=True,
    )

    # 4. Generate small germline VCF for BAF
    vcf_gz = vcf_dir / "mutect2" / "filtered" / "PAT1.filtered.vcf.gz"
    _write_germline_baf_vcf(bgzip, tabix, vcf_gz, seqs, targets)

    # 5. Sample sheet and config (PureCN disabled)
    samples_tsv = tmp_path / "samples.tsv"
    pd.DataFrame(
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
                "sample": "NORM2",
                "tumor_bam": "N2",
                "normal_bam": ".",
                "analysis_type": "germline",
                "in_pon": "yes",
                "sex": "female",
            },
            {
                "sample": "NORM3",
                "tumor_bam": "N3",
                "normal_bam": ".",
                "analysis_type": "germline",
                "in_pon": "yes",
                "sex": "female",
            },
        ]
    ).to_csv(samples_tsv, sep="\t", index=False)

    config_dict: dict[str, Any] = {
        "ref": {
            "genome": str(fasta),
            "build": "GRCh38",
            "targets_bed": str(targets_bed),
            "baits_bed": "",
            "access_bed": "",
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
            "file_extension": ext,
        },
        "cnvkit": {
            "reference_cnn": "",
            "autobin": False,
            "target_avg_size": 200,
            "antitarget_avg_size": 4000,
            "antitarget_min_size": 1000,
            "access_min_gap_size": 1000,
            "pon_bivar_threshold": 0.5,
            "male_reference": False,
            "segment_method": "cbs",
            "segment_threshold": 0.001,
            "drop_low_coverage": True,
            "call_method": "threshold",
            "call_filter": "ci",
            "min_variant_depth": 20,
            "genes_of_interest": ["EGFR", "TP53"],
        },
        "purecn": {
            "enabled": False,
            "mode": "cnvkit",
            "genome": "hg38",
            "assay": "exome",
            "fun_segmentation": "PSCBS",
            "normal_panel_vcf": "",
            "normaldb": "",
            "mapping_bias": "",
            "snp_blacklist": "",
            "seed": 123,
            "postoptimize": True,
            "extra": "",
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config_dict, sort_keys=False), encoding="utf-8")

    # 6. Run real Snakemake workflow (pooled normal reference mode)
    res = subprocess.run(
        [
            snakemake,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(config_path),
            "--directory",
            str(tmp_path),
            "--workflow-profile",
            "none",
            "--cores",
            "2",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=240,
    )
    assert res.returncode == 0, f"Snakemake e2e run failed:\n{res.stdout}\n{res.stderr}"

    # 7. Assert per-sample .cnr/.cns/.call.cns, exports, plots, and QC exist
    cnr_file = out_dir / "cnvkit" / "runs" / "PAT1.cnr"
    cns_file = out_dir / "cnvkit" / "runs" / "PAT1.cns"
    segmetrics_file = out_dir / "cnvkit" / "runs" / "PAT1.segmetrics.cns"
    call_cns_file = out_dir / "cnvkit" / "calls" / "PAT1.call.cns"
    bintest_file = out_dir / "cnvkit" / "calls" / "PAT1.bintest.cns"
    seg_export = out_dir / "cnvkit" / "exports" / "PAT1.seg"
    bed_export = out_dir / "cnvkit" / "exports" / "PAT1.cnv.bed"
    vcf_export = out_dir / "cnvkit" / "exports" / "PAT1.cnv.vcf.gz"
    tbi_export = out_dir / "cnvkit" / "exports" / "PAT1.cnv.vcf.gz.tbi"
    summary_tsv = out_dir / "qc" / "pipeline_summary.tsv"
    normal_qc_tsv = out_dir / "qc" / "normal_qc_summary.tsv"
    clean_normals = out_dir / "cnvkit" / "pon_qc" / "clean_normals.list"

    for path in (
        cnr_file,
        cns_file,
        segmetrics_file,
        call_cns_file,
        bintest_file,
        seg_export,
        bed_export,
        vcf_export,
        tbi_export,
        summary_tsv,
        normal_qc_tsv,
        clean_normals,
        out_dir / "cnvkit" / "reference" / "pooled_reference.cnn",
        out_dir / "cnvkit" / "plots" / "scatter_genome" / "PAT1.genome.pdf",
        out_dir / "cnvkit" / "plots" / "diagram" / "PAT1.diagram.pdf",
        out_dir / "cnvkit" / "plots" / "scatter_gene" / "PAT1" / "EGFR.pdf",
        out_dir / "cnvkit" / "plots" / "scatter_gene" / "PAT1" / "TP53.pdf",
        out_dir / "cnvkit" / "plots" / "heatmap_cohort.pdf",
    ):
        assert path.is_file() and path.stat().st_size > 0, f"Missing or empty output: {path}"

    # Verify all 3 normals passed PoN QC (3 target + 3 antitarget = 6 .cnn paths)
    clean_lines = [ln for ln in clean_normals.read_text(encoding="utf-8").splitlines() if ln]
    assert len(clean_lines) == 6, f"Expected 6 clean normal .cnn files, got: {clean_lines}"

    # 8. Assert simulated gain and loss regions are called with the correct sign and CN
    def _assert_gain_and_loss_calls(call_path: Path, cnr_path: Path) -> None:
        cnr_df = pd.read_csv(cnr_path, sep="\t")
        gain_bins = cnr_df[
            (cnr_df["chromosome"] == GAIN_REGION[0])
            & (cnr_df["start"] >= GAIN_REGION[1])
            & (cnr_df["end"] <= GAIN_REGION[2])
        ]
        loss_bins = cnr_df[
            (cnr_df["chromosome"] == LOSS_REGION[0])
            & (cnr_df["start"] >= LOSS_REGION[1])
            & (cnr_df["end"] <= LOSS_REGION[2])
        ]
        assert gain_bins["log2"].median() > 0.8, f"Expected positive gain log2: {gain_bins}"
        assert loss_bins["log2"].median() < -0.8, f"Expected negative loss log2: {loss_bins}"

        call_df = pd.read_csv(call_path, sep="\t")
        gain_segs = call_df[call_df["chromosome"] == GAIN_REGION[0]]
        loss_segs = call_df[call_df["chromosome"] == LOSS_REGION[0]]
        neut_segs = call_df[call_df["chromosome"].isin(["chr3", "chrX"])]
        assert not gain_segs.empty and (gain_segs["log2"] > 0.5).any()
        assert (gain_segs["cn"] > 2).any(), f"Expected gain cn > 2:\n{gain_segs}"
        assert not loss_segs.empty and (loss_segs["log2"] < -0.5).any()
        assert (loss_segs["cn"] < 2).any(), f"Expected loss cn < 2:\n{loss_segs}"
        assert not neut_segs.empty and (neut_segs["cn"] == 2).all(), (
            f"Expected neutral cn == 2:\n{neut_segs}"
        )

    _assert_gain_and_loss_calls(call_cns_file, cnr_file)

    vcf_view = subprocess.run(
        [bcftools, "view", "-H", str(vcf_export)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "SVTYPE=DUP" in vcf_view.stdout, f"Expected DUP in exported VCF:\n{vcf_view.stdout}"
    assert "SVTYPE=DEL" in vcf_view.stdout, f"Expected DEL in exported VCF:\n{vcf_view.stdout}"

    bed_text = bed_export.read_text(encoding="utf-8")
    assert "chr1\t" in bed_text and "chr2\t" in bed_text

    # 9. Verify the flat-reference fallback path works when no clean normals are available
    flat_out_dir = tmp_path / "results_flat"
    flat_cfg = copy.deepcopy(config_dict)
    flat_cfg["paths"]["output_folder"] = str(flat_out_dir)
    flat_cfg["cnvkit"]["pon_bivar_threshold"] = 0.0
    flat_cfg_path = tmp_path / "config_flat.yaml"
    flat_cfg_path.write_text(yaml.safe_dump(flat_cfg, sort_keys=False), encoding="utf-8")

    # Reuse already-computed per-BAM coverages
    (flat_out_dir / "cnvkit").mkdir(parents=True)
    (flat_out_dir / "cnvkit" / "coverage").symlink_to(out_dir / "cnvkit" / "coverage")

    flat_call_cns = flat_out_dir / "cnvkit" / "calls" / "PAT1.call.cns"
    flat_res = subprocess.run(
        [
            snakemake,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(flat_cfg_path),
            "--directory",
            str(tmp_path),
            "--workflow-profile",
            "none",
            "--cores",
            "2",
            str(flat_call_cns),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert flat_res.returncode == 0, (
        f"Flat-reference fallback run failed:\n{flat_res.stdout}\n{flat_res.stderr}"
    )

    flat_clean = flat_out_dir / "cnvkit" / "pon_qc" / "clean_normals.list"
    assert flat_clean.is_file() and flat_clean.read_text(encoding="utf-8") == ""
    flat_ref_log = (flat_out_dir / "logs" / "cnvkit_reference_pooled" / "reference.log").read_text(
        encoding="utf-8"
    )
    assert "WARNING: No clean normals available; building flat reference." in flat_ref_log
    assert (flat_out_dir / "cnvkit" / "reference" / "pooled_reference.cnn").is_file()
    _assert_gain_and_loss_calls(
        flat_call_cns,
        flat_out_dir / "cnvkit" / "runs" / "PAT1.cnr",
    )
