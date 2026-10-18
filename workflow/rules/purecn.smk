"""PureCN tumor purity and ploidy estimation rules (cnvkit and native modes)."""

if PURECN_ENABLED:

    rule purecn_mapping_bias:
        """Build PureCN mapping bias RDS from a normal panel VCF (cnvkit mode)."""
        input:
            normal_panel=PURECN_NORMAL_PANEL_VCF if PURECN_NORMAL_PANEL_VCF else [],
        output:
            bias=os.path.join(PURECN_REF_DIR, f"mapping_bias_{PURECN_ASSAY}_{PURECN_GENOME}.rds"),
        params:
            outdir=PURECN_REF_DIR,
            assay=PURECN_ASSAY,
            genome=PURECN_GENOME,
            extra=PURECN_NORMALDB_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wc, attempt: 8000 * attempt,
            runtime=lambda wc, attempt: min(240 * attempt, 480),
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_mapping_bias/mapping_bias.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/purecn_mapping_bias/mapping_bias.tsv")
        shell:
            r"""
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')
            Rscript "$PURECN/NormalDB.R" \
                --outdir {params.outdir} \
                --normal-panel {input.normal_panel} \
                --assay {params.assay} \
                --genome {params.genome} \
                {params.extra} \
                --force &> {log}
            """

    rule purecn_intervals:
        """Generate GC- and mappability-calibrated interval file for PureCN (native mode)."""
        input:
            bed=lambda wc: (
                os.path.join(CNVKIT_REF_DIR, "targets.clean.bed")
                if REF_BLACKLIST_BED and BAITS_BED == TARGETS_BED
                else BAITS_BED
            ),
            fasta=REF,
            mappability=[REF_MAPPABILITY_BW] if REF_MAPPABILITY_BW else [],
            reptiming=[REF_REPTIMING_BW] if REF_REPTIMING_BW else [],
        output:
            intervals=os.path.join(PURECN_REF_DIR, "intervals.txt"),
        params:
            genome=PURECN_GENOME,
            mappability_arg=(f"--mappability {REF_MAPPABILITY_BW}" if REF_MAPPABILITY_BW else ""),
            reptiming_arg=(f"--reptiming {REF_REPTIMING_BW}" if REF_REPTIMING_BW else ""),
            extra=PURECN_INTERVALFILE_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wc, attempt: 8000 * attempt,
            runtime=lambda wc, attempt: min(120 * attempt, 480),
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_intervals/intervals.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/purecn_intervals/intervals.tsv")
        shell:
            r"""
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')
            Rscript "$PURECN/IntervalFile.R" \
                --infile {input.bed} \
                --fasta {input.fasta} \
                --outfile {output.intervals} \
                --offtarget \
                --genome {params.genome} \
                {params.mappability_arg} \
                {params.reptiming_arg} \
                {params.extra} &> {log}
            """

    rule purecn_coverage:
        """Compute GC-normalized PureCN coverage per unique BAM (native mode)."""
        input:
            bam=get_bam_by_basename,
            intervals=os.path.join(PURECN_REF_DIR, "intervals.txt"),
        output:
            coverage=os.path.join(PURECN_COV_DIR, "{bam_basename}_coverage_loess.txt.gz"),
        params:
            outdir=PURECN_COV_DIR,
            extra=PURECN_COVERAGE_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wc, attempt: 8000 * attempt,
            runtime=lambda wc, attempt: min(240 * attempt, 480),
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_coverage/{bam_basename}.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/purecn_coverage/{bam_basename}.tsv")
        shell:
            r"""
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')
            Rscript "$PURECN/Coverage.R" \
                --bam {input.bam} \
                --intervals {input.intervals} \
                --outdir {params.outdir} \
                {params.extra} &> {log}
            """

    rule purecn_normaldb:
        """Build PureCN normalDB.rds and mapping_bias.rds from normal BAM coverages (native mode)."""
        input:
            coverages=lambda wc: [
                os.path.join(PURECN_COV_DIR, f"{bn}_coverage_loess.txt.gz")
                for bn in get_normal_bams()
            ],
            normal_panel=[PURECN_NORMAL_PANEL_VCF] if PURECN_NORMAL_PANEL_VCF else [],
        output:
            normaldb=os.path.join(PURECN_REF_DIR, "normalDB.rds"),
            bias=os.path.join(PURECN_REF_DIR, "mapping_bias.rds"),
        params:
            outdir=PURECN_REF_DIR,
            genome=PURECN_GENOME,
            normal_panel_arg=(
                f"--normal-panel {PURECN_NORMAL_PANEL_VCF}" if PURECN_NORMAL_PANEL_VCF else ""
            ),
            extra=PURECN_NORMALDB_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wc, attempt: 16000 * attempt,
            runtime=lambda wc, attempt: min(480 * attempt, 480),
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_normaldb/normaldb.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/purecn_normaldb/normaldb.tsv")
        shell:
            r"""
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')
            Rscript "$PURECN/NormalDB.R" \
                --coveragefiles {input.coverages} \
                --outdir {params.outdir} \
                --genome {params.genome} \
                {params.normal_panel_arg} \
                {params.extra} \
                --force &> {log}
            if [ ! -f "{output.bias}" ]; then
                FOUND_BIAS=$(ls {params.outdir}/mapping_bias*.rds 2>/dev/null | head -n 1 || true)
                if [ -n "$FOUND_BIAS" ]; then
                    cp "$FOUND_BIAS" "{output.bias}"
                else
                    touch "{output.bias}"
                fi
            fi
            """

    def _purecn_run_tumor_input(wildcards):
        if PURECN_MODE == "native":
            bn = get_unit_tumor_bam_basename(wildcards.unit)
            return os.path.join(PURECN_COV_DIR, f"{bn}_coverage_loess.txt.gz")
        return os.path.join(CNVKIT_RUNS_DIR, f"{wildcards.unit}.cnr")

    def _purecn_run_seg_input(wildcards):
        if PURECN_MODE == "cnvkit":
            return [os.path.join(CNVKIT_EXPORTS_DIR, f"{wildcards.unit}.seg")]
        return []

    def _purecn_run_intervals_input(wildcards):
        if PURECN_MODE == "native":
            return [os.path.join(PURECN_REF_DIR, "intervals.txt")]
        return []

    def _purecn_run_normaldb_input(wildcards):
        if PURECN_MODE == "native":
            return [get_purecn_normaldb_path()]
        return []

    def _purecn_run_mapping_bias_input(wildcards):
        mb = get_purecn_mapping_bias_path()
        return [mb] if mb else []

    def _purecn_run_stats_input(wildcards):
        st = get_unit_vcf_stats_path(wildcards.unit)
        return [st] if st else []

    def _purecn_run_mode_args(wildcards, input):
        parts = [f"--tumor {input.tumor}"]
        if PURECN_MODE == "native":
            parts.append(f"--intervals {input.intervals[0]}")
            parts.append(f"--normaldb {input.normaldb[0]}")
        else:
            parts.append(f"--seg-file {input.seg[0]}")
            parts.append(f"--fun-segmentation {PURECN_FUN_SEGMENTATION}")
        if input.mapping_bias:
            parts.append(f"--mappingbiasfile {input.mapping_bias[0]}")
        return " ".join(parts)

    rule purecn_run:
        """Run fault-tolerant PureCN purity and ploidy estimation for a tumor unit."""
        input:
            tumor=_purecn_run_tumor_input,
            seg=_purecn_run_seg_input,
            intervals=_purecn_run_intervals_input,
            normaldb=_purecn_run_normaldb_input,
            mapping_bias=_purecn_run_mapping_bias_input,
            vcf=lambda wc: get_unit_vcf_path(wc.unit),
            stats=_purecn_run_stats_input,
            snp_blacklist=[PURECN_SNP_BLACKLIST] if PURECN_SNP_BLACKLIST else [],
        output:
            csv=os.path.join(PURECN_RUNS_DIR, "{unit}", "{unit}.csv"),
            rds=os.path.join(PURECN_RUNS_DIR, "{unit}", "{unit}.rds"),
            status=os.path.join(PURECN_RUNS_DIR, "{unit}", "{unit}.status"),
        params:
            outdir=lambda wc: os.path.join(PURECN_RUNS_DIR, wc.unit),
            sampleid=lambda wc: CNV_UNITS[wc.unit]["tumor_vcf_id"],
            mode_args=_purecn_run_mode_args,
            stats_arg=lambda wc, input: (f"--statsfile {input.stats[0]}" if input.stats else ""),
            snp_blacklist_arg=(
                f"--snpblacklist {PURECN_SNP_BLACKLIST}" if PURECN_SNP_BLACKLIST else ""
            ),
            genome=PURECN_GENOME,
            seed=PURECN_SEED,
            postoptimize_arg="--postoptimize" if PURECN_POSTOPTIMIZE else "",
            extra=PURECN_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wc, attempt: 16000 * attempt,
            runtime=lambda wc, attempt: min(480 * attempt, 480),
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_run/{unit}.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/purecn_run/{unit}.tsv")
        shell:
            r"""
            mkdir -p {params.outdir}
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')
            if Rscript "$PURECN/PureCN.R" \
                --out {params.outdir}/{wildcards.unit} \
                --sampleid {params.sampleid} \
                {params.mode_args} \
                --vcf {input.vcf} \
                {params.stats_arg} \
                --genome {params.genome} \
                --force \
                --seed {params.seed} \
                {params.postoptimize_arg} \
                {params.snp_blacklist_arg} \
                {params.extra} \
                > {log} 2>&1; then
                echo "succeeded" > {output.status}
                touch {output.csv} {output.rds}
            else
                echo "WARNING: PureCN.R failed to converge for {wildcards.unit}; falling back." >> {log}
                echo "failed" > {output.status}
                : > {output.csv}
                : > {output.rds}
            fi
            """


rule consolidate_purity:
    """Consolidate tumor purity and ploidy from PureCN or sample sheet metadata."""
    localrule: True
    input:
        csv=lambda wc: (
            [os.path.join(PURECN_RUNS_DIR, wc.unit, f"{wc.unit}.csv")]
            if wc.unit in get_purecn_units()
            else []
        ),
        status=lambda wc: (
            [os.path.join(PURECN_RUNS_DIR, wc.unit, f"{wc.unit}.status")]
            if wc.unit in get_purecn_units()
            else []
        ),
    output:
        tsv=os.path.join(PURITY_DIR, "{unit}.purity_ploidy.tsv"),
    params:
        sample=lambda wc: CNV_UNITS[wc.unit]["sample"],
        purecn_csv=lambda wc, input: input.csv[0] if input.csv else "",
        purecn_status=lambda wc, input: input.status[0] if input.status else "",
        sample_purity=lambda wc: (
            "" if CNV_UNITS[wc.unit]["purity"] is None else CNV_UNITS[wc.unit]["purity"]
        ),
        sample_ploidy=lambda wc: (
            "" if CNV_UNITS[wc.unit]["ploidy"] is None else CNV_UNITS[wc.unit]["ploidy"]
        ),
        purecn_attempted=lambda wc: wc.unit in get_purecn_units(),
        rules_dir=RULES_DIR,
    log:
        os.path.join(LOG_DIR, "consolidate_purity/{unit}.log"),
    shell:
        r"""
        python3 -c '
import sys
sys.path.insert(0, "{params.rules_dir}")
from helpers import consolidate_purity_ploidy
consolidate_purity_ploidy(
    unit="{wildcards.unit}",
    sample="{params.sample}",
    output_file="{output.tsv}",
    purecn_csv="{params.purecn_csv}",
    purecn_status_file="{params.purecn_status}",
    sample_purity="{params.sample_purity}",
    sample_ploidy="{params.sample_ploidy}",
    purecn_attempted=bool({params.purecn_attempted}),
)
' &> {log}
        """
