import os
import pandas as pd

# ----------------------------------------------------------------------------------- #
#                                 CNVkit-PureCN Pipeline                              #
# ----------------------------------------------------------------------------------- #
# This Snakemake workflow automates copy number variation (CNV) analysis from
# tumor exome data, implementing best practices for noisy samples.
# See README.md for a full description.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Configuration and Setup
# ----------------------------------------------------------------------------------- #
configfile: "config.yaml"

# Load samplesheet into a pandas DataFrame, keyed by sample_id for easy lookup
try:
    SAMPLES = pd.read_csv(config["samplesheet"], sep="\t", dtype=str).set_index("sample_id", drop=False)
except FileNotFoundError:
    raise FileNotFoundError(f"Sample sheet not found at {config['samplesheet']}. Please create it.")

# Helper functions to get lists of sample IDs based on analysis type
def get_tumor_samples():
    """Returns a list of all tumor sample IDs."""
    return SAMPLES[SAMPLES.analysis_type.isin(["TvsN", "To"])].index.tolist()

def get_normal_sample_ids():
    """Returns a list of all normal sample IDs."""
    return SAMPLES[SAMPLES.analysis_type == "Normal"].index.tolist()

def get_purecn_candidates():
    """Returns sample_ids eligible for PureCN (TvsN with a valid VCF path)."""
    return SAMPLES[
        (SAMPLES.analysis_type == "TvsN") & (SAMPLES.vcf.notna()) & (SAMPLES.vcf != "")
    ].index.tolist()

def get_vcf_for_call(wildcards):
    """Returns VCF path for cnvkit_call if available, empty list otherwise."""
    sample_vcf = SAMPLES.loc[wildcards.sample_id, "vcf"]
    if pd.notna(sample_vcf) and sample_vcf != "":
        return [sample_vcf]
    return []

# ----------------------------------------------------------------------------------- #
# Rule 'all': Defines the final targets of the workflow
# ----------------------------------------------------------------------------------- #
rule all:
    input:
        expand("results/07_final_vcfs/{sample_id}.cnv.vcf.gz", sample_id=get_tumor_samples())

# ----------------------------------------------------------------------------------- #
# MODULE 1: PURITY ESTIMATION
# ----------------------------------------------------------------------------------- #
rule create_purecn_mapping_bias:
    output:
        mapping_bias_db=f"results/01_purecn_setup/mapping_bias_{config['purecn_assay_name']}_{config['purecn_genome']}.rds"
    params:
        assay=config["purecn_assay_name"],
        genome=config["purecn_genome"],
        normal_panel_vcf=config["purecn_normal_panel_vcf"]
    log:
        "logs/purecn/create_mapping_bias.log"
    conda:
        f"{config['conda_env_dir']}/purecn.yaml"
    shell:
        """
        Rscript $(Rscript -e "cat(system.file('extdata', 'NormalDB.R', package='PureCN'))") \\
            --out-dir results/01_purecn_setup \\
            --normal-panel {params.normal_panel_vcf} \\
            --assay {params.assay} \\
            --genome {params.genome} \\
            --force &> {log}
        """

rule run_purecn:
    input:
        tumor_bam=lambda w: SAMPLES.loc[w.sample_id, "tumor_bam"],
        vcf=lambda w: SAMPLES.loc[w.sample_id, "vcf"],
        mapping_bias_db=rules.create_purecn_mapping_bias.output.mapping_bias_db
    output:
        rds="results/02_purecn_runs/{sample_id}.rds",
        csv="results/02_purecn_runs/{sample_id}.csv"
    params:
        genome=config["purecn_genome"],
        out_dir="results/02_purecn_runs"
    log:
        "logs/purecn/run_purecn.{sample_id}.log"
    conda:
        f"{config['conda_env_dir']}/purecn.yaml"
    threads: config["default_threads"]
    resources:
        mem_mb=config["default_mem_mb"]
    shell:
        """
        # Create sample-specific temp directory
        mkdir -p {params.out_dir}/temp_{wildcards.sample_id}
        
        # PureCN requires a .cnr file. We generate a temporary one using a flat reference.
        # This provides the necessary coverage info without needing the final PoN.
        cnvkit.py coverage {input.tumor_bam} {config[targets_bed]} -p {threads} -o {params.out_dir}/temp_{wildcards.sample_id}/temp.target.cnn &>> {log}
        cnvkit.py coverage {input.tumor_bam} {config[access_bed]} -p {threads} -o {params.out_dir}/temp_{wildcards.sample_id}/temp.antitarget.cnn &>> {log}
        cnvkit.py reference -f {config[reference_genome]} -t {config[targets_bed]} -a {config[access_bed]} -o {params.out_dir}/temp_{wildcards.sample_id}/temp.ref.cnn &>> {log}
        cnvkit.py fix {params.out_dir}/temp_{wildcards.sample_id}/temp.target.cnn {params.out_dir}/temp_{wildcards.sample_id}/temp.antitarget.cnn {params.out_dir}/temp_{wildcards.sample_id}/temp.ref.cnn -o {params.out_dir}/temp_{wildcards.sample_id}/temp.cnr &>> {log}

        # Now run PureCN with the temporary CNR file
        Rscript $(Rscript -e "cat(system.file('extdata', 'PureCN.R', package='PureCN'))") \\
            --out {params.out_dir}/{wildcards.sample_id} \\
            --sampleid {wildcards.sample_id} \\
            --tumor {params.out_dir}/temp_{wildcards.sample_id}/temp.cnr \\
            --vcf {input.vcf} \\
            --mapping-bias-file {input.mapping_bias_db} \\
            --genome {params.genome} \\
            --post-optimize --force --seed 123 &>> {log}
        
        # Move outputs to final locations
        mv {params.out_dir}/{wildcards.sample_id}.csv {output.csv}
        mv {params.out_dir}/{wildcards.sample_id}.rds {output.rds}
        
        # Clean up temporary files
        rm -rf {params.out_dir}/temp_{wildcards.sample_id}
        """

rule consolidate_purity:
    input:
        purecn_csv=lambda w: f"results/02_purecn_runs/{w.sample_id}.csv" if w.sample_id in get_purecn_candidates() else [],
    output:
        purity_file="results/03_purity_values/{sample_id}.purity.txt"
    params:
        sample_id="{sample_id}",
        fallback_purity=lambda w: SAMPLES.loc[w.sample_id, "fallback_purity"],
        # Pass a string 'None' if the optional input file does not exist
        purecn_csv_param=lambda w, input: input.purecn_csv[0] if input.purecn_csv else "None"
    log:
        "logs/purity/consolidate_purity.{sample_id}.log"
    shell:
        """
        python scripts/snakemake/helpers/consolidate_purity.py \\
            --purecn-csv {params.purecn_csv_param} \\
            --fallback-purity {params.fallback_purity} \\
            --output {output.purity_file} &> {log}
        """

# ----------------------------------------------------------------------------------- #
# MODULE 2: PANEL OF NORMALS (PoN) GENERATION
# ----------------------------------------------------------------------------------- #
rule cnvkit_coverage_normals:
    input:
        bam=lambda w: SAMPLES.loc[w.normal_id, "normal_bam"]
    output:
        target="results/04_pon_creation/coverage/{normal_id}.targetcoverage.cnn",
        antitarget=temp("results/04_pon_creation/coverage/{normal_id}.antitargetcoverage.cnn")
    log:
        "logs/pon/cnvkit_coverage.{normal_id}.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    threads: config["default_threads"]
    shell:
        """
        cnvkit.py coverage {input.bam} {config[targets_bed]} -p {threads} -o {output.target} &> {log}
        # Antitargets are intermediate and can be temp
        cnvkit.py coverage {input.bam} {config[access_bed]} -p {threads} -o {output.antitarget} &>> {log}
        """

rule cnvkit_metrics_normals:
    input:
        # Collect all target coverage files for normal samples
        expand(rules.cnvkit_coverage_normals.output.target, normal_id=get_normal_sample_ids())
    output:
        metrics_file="results/04_pon_creation/qc/all_normals.metrics.txt"
    log:
        "logs/pon/cnvkit_metrics_normals.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        """
        # The 'sample' column will contain the full path to the input .cnn files
        cnvkit.py metrics {input} > {output.metrics_file} 2> {log}
        """

rule identify_clean_normals:
    input:
        metrics=rules.cnvkit_metrics_normals.output.metrics_file
    output:
        clean_list="results/04_pon_creation/qc/clean_normals.list"
    params:
        threshold=config["pon_qc_metric_threshold"]
    log:
        "logs/pon/identify_clean_normals.log"
    shell:
        """
        python scripts/snakemake/helpers/identify_clean_normals.py -i {input.metrics} -o {output.clean_list} -t {params.threshold} &> {log}
        """

rule cnvkit_reference_pooled:
    input:
        normal_list=rules.identify_clean_normals.output.clean_list,
        # This ensures all normal coverages are generated before this rule tries to use them
        normal_coverages=expand(rules.cnvkit_coverage_normals.output.target, normal_id=get_normal_sample_ids())
    output:
        pooled_ref="results/04_pon_creation/pooled_reference.cnn"
    log:
        "logs/pon/cnvkit_reference_pooled.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        # Create the glob patterns from the list of clean normals
        # e.g., /path/to/normal1.targetcoverage.cnn -> /path/to/normal1.*targetcoverage.cnn
        # This will match both target and antitarget files for the clean normals.
        """
        if [ -s {input.normal_list} ]; then
            GLOB_PATTERN=$(cat {input.normal_list} | sed 's/targetcoverage/\\*targetcoverage/' | tr '\\n' ' ')
            cnvkit.py reference $GLOB_PATTERN \\
                -f {config[reference_genome]} \\
                -o {output.pooled_ref} &> {log}
        else
            echo "WARNING: No clean normals found. Creating a flat reference instead." > {log}
            cnvkit.py reference -f {config[reference_genome]} \\
                -t {config[targets_bed]} \\
                -a {config[access_bed]} \\
                -o {output.pooled_ref} &>> {log}
        fi
        """

# ----------------------------------------------------------------------------------- #
# MODULE 3: CNV CALLING IN TUMOR SAMPLES
# ----------------------------------------------------------------------------------- #
rule cnvkit_batch_tumor:
    input:
        tumor_bam=lambda w: SAMPLES.loc[w.sample_id, "tumor_bam"],
        pooled_ref=rules.cnvkit_reference_pooled.output.pooled_ref
    output:
        cnr="results/05_cnvkit_runs/{sample_id}.cnr",
        cns="results/05_cnvkit_runs/{sample_id}.cns"
    log:
        "logs/cnvkit/cnvkit_batch.{sample_id}.log"
    params:
        output_dir="results/05_cnvkit_runs/"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    threads: config["default_threads"]
    shell:
        """
        cnvkit.py batch {input.tumor_bam} \\
            -r {input.pooled_ref} \\
            -p {threads} \\
            --exclude {config[plasmid_blacklist]} \\
            --segment-threshold {config[segment_threshold]} \\
            -d {params.output_dir} &> {log}
        """

rule cnvkit_call:
    input:
        cns=rules.cnvkit_batch_tumor.output.cns,
        purity_file=rules.consolidate_purity.output.purity_file,
        vcf=get_vcf_for_call, # Helper function handles optional input
    output:
        call_cns="results/06_final_calls/{sample_id}.call.cns"
    log:
        "logs/cnvkit/cnvkit_call.{sample_id}.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        """
        # Conditionally add the -v flag only if the input.vcf list is not empty
        VCF_PARAM=""
        if [ -n "{input.vcf}" ]; then
            VCF_PARAM="-v {input.vcf}"
        fi
        
        cnvkit.py call {input.cns} \\
            $VCF_PARAM \\
            --purity $(cat {input.purity_file}) \\
            -m clonal \\
            -o {output.call_cns} &> {log}
        """

# ----------------------------------------------------------------------------------- #
# MODULE 4: FINAL EXPORT
# ----------------------------------------------------------------------------------- #
rule cnvkit_export_vcf:
    input:
        call_cns=rules.cnvkit_call.output.call_cns
    output:
        vcf="results/07_final_vcfs/{sample_id}.cnv.vcf.gz",
        tbi="results/07_final_vcfs/{sample_id}.cnv.vcf.gz.tbi"
    params:
        sample_id="{sample_id}"
    log:
        "logs/cnvkit/cnvkit_export.{sample_id}.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        """
        cnvkit.py export vcf {input.call_cns} -i {params.sample_id} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
