# ----------------------------------------------------------------------------------- #
# MODULE 2: PANEL OF NORMALS (PoN) GENERATION
# ----------------------------------------------------------------------------------- #

ruleorder: cnvkit_reference_pooled > cnvkit_batch_tumor

rule cnvkit_coverage_normals:
    input:
        bam=lambda w: get_normal_bams_for_pon()[w.normal_id]
    output:
        target=f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.targetcoverage.cnn",
        antitarget=temp(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.antitargetcoverage.cnn")
    log:
        f"{config['dirs']['logs']}/cnvkit_coverage_normals/{{normal_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
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
        config["conda_envs"]["cnvkit"]
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
        normal_coverages=expand(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.targetcoverage.cnn", normal_id=get_normal_sample_ids())
    output:
        pooled_ref=f"{config['dirs']['pon_creation']}/pooled_reference.cnn"
    params:
        # Conditionally add annotate flag if refFlat is provided in config
        annotate_flag=lambda wildcards: f"--annotate {config['annotate_refFlat']}" if "annotate_refFlat" in config and config["annotate_refFlat"] else ""
    log:
        f"{config['dirs']['logs']}/cnvkit_reference_pooled/log.txt"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        """
        # Note: The 'cnvkit.py batch' command uses '--targets' and '--annotate' to create an
        # intermediate, annotated target file. The 'cnvkit.py reference' command, however,
        # does not have a direct '--annotate' flag. Annotation is part of the 'target'
        # command. To handle this correctly in a modular pipeline, we should create an
        # annotated target file first.
        #
        # CORRECTION: The best practice is to annotate the reference file itself. The
        # `cnvkit reference` command DOES have an --annotate flag when building from
        # coverage files. My previous thinking was incorrect. The logic below is what
        # was intended.
        # Let's adjust the logic slightly for robustness. `cnvkit reference` can be picky.
        # We will build the command string carefully.

        if [ -s {input.normal_list} ]; then
            # Build reference from a panel of clean normals
            echo "Building pooled reference from {input.normal_list}" > {log}
            GLOB_PATTERN=$(cat {input.normal_list} | sed 's/\\.targetcoverage\\.cnn/\\.\\*targetcoverage\\.cnn/' | tr '\\n' ' ')
            cnvkit.py reference $GLOB_PATTERN \\
                -f {config[reference_genome]} \\
                {params.annotate_flag} \\
                -o {output.pooled_ref} &>> {log}
        else
            # Fallback to a flat reference if no clean normals are available
            echo 'WARNING: No clean normals found. Creating a flat reference instead.' > {log}
            cnvkit.py reference \\
                -f {config[reference_genome]} \\
                -t {config[targets_bed]} \\
                -a {config[access_bed]} \\
                {params.annotate_flag} \\
                -o {output.pooled_ref} &>> {log}
        fi
        """
