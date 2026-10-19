"""QC, visualization (scatter, diagram, heatmap), and cohort summary report rules."""


rule cnvkit_metrics_tumor:
    """Compute sample-level dispersion/noise metrics (stdev, mad, iqr, bivar) per tumor unit."""
    input:
        cnr=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cnr"),
        cns=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cns"),
    output:
        metrics=os.path.join(QC_DIR, "tumor_metrics", "{unit}.metrics.tsv"),
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(60 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_metrics_tumor/{unit}.log"),
    shell:
        r"""
        cnvkit.py metrics {input.cnr} -s {input.cns} -o {output.metrics} &> {log}
        """


rule cnvkit_scatter_genome:
    """Generate genome-wide log2 copy-ratio and BAF scatter plot per tumor unit."""
    input:
        cnr=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cnr"),
        cns=os.path.join(CNVKIT_CALLS_DIR, "{unit}.call.cns"),
        vcf=get_unit_vcf_input,
    output:
        pdf=os.path.join(CNVKIT_PLOTS_DIR, "scatter_genome", "{unit}.genome.pdf"),
    params:
        vcf_args=get_unit_vcf_args,
        extra=CNVKIT_SCATTER_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 6000 * attempt,
        runtime=lambda wc, attempt: min(120 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_scatter_genome/{unit}.log"),
    shell:
        r"""
        cnvkit.py scatter {input.cnr} \
            -s {input.cns} \
            {params.vcf_args} \
            {params.extra} \
            -o {output.pdf} &> {log}
        """


rule cnvkit_scatter_gene:
    """Generate gene-focused scatter plot per tumor unit and gene of interest."""
    input:
        cnr=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cnr"),
        cns=os.path.join(CNVKIT_CALLS_DIR, "{unit}.call.cns"),
        vcf=get_unit_vcf_input,
    output:
        pdf=os.path.join(CNVKIT_PLOTS_DIR, "scatter_gene", "{unit}", "{gene}.pdf"),
    params:
        vcf_args=get_unit_vcf_args,
        extra=CNVKIT_SCATTER_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 6000 * attempt,
        runtime=lambda wc, attempt: min(60 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_scatter_gene/{unit}.{gene}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.pdf})
        cnvkit.py scatter {input.cnr} \
            -s {input.cns} \
            {params.vcf_args} \
            -g {wildcards.gene} \
            {params.extra} \
            -o {output.pdf} &> {log}
        """


rule cnvkit_diagram:
    """Generate chromosome ideogram copy-number diagram with sex-chromosome handling."""
    input:
        cnr=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cnr"),
        cns=os.path.join(CNVKIT_CALLS_DIR, "{unit}.call.cns"),
        sex_txt=os.path.join(CNVKIT_CALLS_DIR, "{unit}.sex.txt"),
    output:
        pdf=os.path.join(CNVKIT_PLOTS_DIR, "diagram", "{unit}.diagram.pdf"),
    params:
        male_ref_arg="-y" if CNVKIT_MALE_REFERENCE else "",
        extra=CNVKIT_DIAGRAM_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 6000 * attempt,
        runtime=lambda wc, attempt: min(60 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_diagram/{unit}.log"),
    shell:
        r"""
        SAMPLE_SEX=$(cat {input.sex_txt})
        cnvkit.py diagram {input.cnr} \
            -s {input.cns} \
            {params.male_ref_arg} \
            -x "$SAMPLE_SEX" \
            {params.extra} \
            -o {output.pdf} &> {log}
        """


rule cnvkit_heatmap_cohort:
    """Generate cohort-wide copy-number heatmap across all tumor units."""
    input:
        cns=lambda wc: expand(
            os.path.join(CNVKIT_CALLS_DIR, "{unit}.call.cns"),
            unit=get_cnv_unit_ids(),
        ),
    output:
        pdf=os.path.join(CNVKIT_PLOTS_DIR, "heatmap_cohort.pdf"),
    params:
        extra=CNVKIT_HEATMAP_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 16000 * attempt,
        runtime=lambda wc, attempt: min(120 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_heatmap_cohort/heatmap.log"),
    shell:
        r"""
        cnvkit.py heatmap {input.cns} -d {params.extra} -o {output.pdf} &> {log}
        """


rule create_summary_report:
    """Aggregate cohort purity, ploidy, chromosomal sex, tumor metrics, and normal PoN QC."""
    localrule: True
    input:
        purity_files=lambda wc: expand(
            os.path.join(PURITY_DIR, "{unit}.purity_ploidy.tsv"),
            unit=get_cnv_unit_ids(),
        ),
        sex_files=lambda wc: expand(
            os.path.join(CNVKIT_CALLS_DIR, "{unit}.sex.txt"),
            unit=get_cnv_unit_ids(),
        ),
        tumor_metrics=lambda wc: expand(
            os.path.join(QC_DIR, "tumor_metrics", "{unit}.metrics.tsv"),
            unit=get_cnv_unit_ids(),
        ),
        normal_metrics=lambda wc: (
            [os.path.join(CNVKIT_PON_QC_DIR, "all_normals.metrics.txt")]
            if get_normal_bams() and not CNVKIT_REFERENCE_CNN
            else []
        ),
        clean_normals=lambda wc: (
            [os.path.join(CNVKIT_PON_QC_DIR, "clean_normals.list")]
            if get_normal_bams() and not CNVKIT_REFERENCE_CNN
            else []
        ),
    output:
        summary=os.path.join(QC_DIR, "pipeline_summary.tsv"),
        normal_qc=os.path.join(QC_DIR, "normal_qc_summary.tsv"),
    params:
        units_meta=CNV_UNITS,
        normal_metrics_path=lambda wc, input: (
            input.normal_metrics[0] if input.normal_metrics else ""
        ),
        clean_normals_path=lambda wc, input: (
            input.clean_normals[0] if input.clean_normals else ""
        ),
        rules_dir=RULES_DIR,
    log:
        os.path.join(LOG_DIR, "create_summary_report/summary.log"),
    shell:
        r"""
        python3 - <<'PY' &> {log}
import sys
sys.path.insert(0, {params.rules_dir!r})
from helpers import build_pipeline_summary
build_pipeline_summary(
    units_meta={params.units_meta!r},
    purity_files=list({input.purity_files!r}),
    sex_files=list({input.sex_files!r}),
    tumor_metrics_files=list({input.tumor_metrics!r}),
    normal_metrics_file={params.normal_metrics_path!r} or None,
    clean_normals_file={params.clean_normals_path!r} or None,
    output_summary_tsv={output.summary!r},
    output_normal_qc_tsv={output.normal_qc!r},
)
PY
        """
