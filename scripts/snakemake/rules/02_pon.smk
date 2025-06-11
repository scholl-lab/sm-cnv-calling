# ----------------------------------------------------------------------------------- #
# MODULE 2: PANEL OF NORMALS (PoN) GENERATION
# ----------------------------------------------------------------------------------- #

ruleorder: cnvkit_reference_pooled > cnvkit_batch_tumor

rule cnvkit_coverage_normals:
    input:
        bam=lambda w: SAMPLES.loc[w.normal_id, "normal_bam"]
    output:
        target=f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.targetcoverage.cnn",
        antitarget=temp(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.antitargetcoverage.cnn")
    log:
        f"{config['dirs']['logs']}/cnvkit_coverage_normals/{{normal_id}}.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    threads: config["default_threads"]
    shell:
        """
        cnvkit.py coverage {input.bam} {config[targets_bed]} -p {threads} -o {output.target} &> {log}
        cnvkit.py coverage {input.bam} {config[access_bed]} -p {threads} -o {output.antitarget} &>> {log}
        """

rule cnvkit_metrics_normals:
    input:
        expand(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.targetcoverage.cnn", normal_id=get_normal_sample_ids())
    output:
        metrics_file=f"{config['dirs']['pon_creation']}/qc/all_normals.metrics.txt"
    log:
        f"{config['dirs']['logs']}/cnvkit_metrics_normals/log.txt"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        "cnvkit.py metrics {input} > {output.metrics_file} 2> {log}"

rule identify_clean_normals:
    input:
        metrics=rules.cnvkit_metrics_normals.output.metrics_file
    output:
        clean_list=f"{config['dirs']['pon_creation']}/qc/clean_normals.list"
    params:
        threshold=config["pon_qc_metric_threshold"]
    log:
        f"{config['dirs']['logs']}/identify_clean_normals/log.txt"
    shell:
        "python scripts/snakemake/helpers/identify_clean_normals.py -i {input.metrics} -o {output.clean_list} -t {params.threshold} &> {log}"

rule cnvkit_reference_pooled:
    input:
        normal_list=rules.identify_clean_normals.output.clean_list,
        normal_coverages=expand(rules.cnvkit_coverage_normals.output.target, normal_id=get_normal_sample_ids())
    output:
        pooled_ref=f"{config['dirs']['pon_creation']}/pooled_reference.cnn"
    log:
        f"{config['dirs']['logs']}/cnvkit_reference_pooled/log.txt"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        """
        if [ -s {input.normal_list} ]; then
            GLOB_PATTERN=$(cat {input.normal_list} | sed 's/targetcoverage/\\*targetcoverage/' | tr '\\n' ' ')
            cnvkit.py reference $GLOB_PATTERN -f {config[reference_genome]} -o {output.pooled_ref} &> {log}
        else
            echo 'WARNING: No clean normals found. Creating a flat reference instead.' > {log}
            cnvkit.py reference -f {config[reference_genome]} -t {config[targets_bed]} -a {config[access_bed]} -o {output.pooled_ref} &>> {log}
        fi
        """
