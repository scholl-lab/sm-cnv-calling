# ----------------------------------------------------------------------------------- #
# MODULE 2: PANEL OF NORMALS (PoN) GENERATION
# ----------------------------------------------------------------------------------- #

ruleorder: cnvkit_reference_pooled > cnvkit_batch_tumor

# Global rule to create antitarget regions (needed by both PoN and PureCN workflows)
rule create_antitarget_bed:
    input:
        targets=config["targets_bed"],
        access=config["access_bed"]
    output:
        antitargets=f"{config['dirs']['pon_creation']}/antitargets.bed"
    log:
        f"{config['dirs']['logs']}/create_antitarget_bed/log.txt"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        "cnvkit.py antitarget {input.targets} -g {input.access} -o {output.antitargets} &> {log}"

rule cnvkit_coverage_normals:
    input:
        bam=lambda w: get_normal_bams_for_pon()[w.normal_id],
        targets=config["targets_bed"],
        antitargets=rules.create_antitarget_bed.output.antitargets
    output:
        target=f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.targetcoverage.cnn",
        antitarget=f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.antitargetcoverage.cnn"
    log:
        f"{config['dirs']['logs']}/cnvkit_coverage_normals/{{normal_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    threads: config["default_threads"]
    shell:
        """
        cnvkit.py coverage {input.bam} {input.targets} -p {threads} -o {output.target} &> {log}
        cnvkit.py coverage {input.bam} {input.antitargets} -p {threads} -o {output.antitarget} &>> {log}
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
        "python {config[helpers_dir]}/identify_clean_normals.py -i {input.metrics} -o {output.clean_list} -t {params.threshold} &> {log}"

rule cnvkit_reference_pooled:
    input:
        normal_list=rules.identify_clean_normals.output.clean_list,
        target_coverages=expand(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.targetcoverage.cnn", normal_id=get_normal_sample_ids()),
        antitarget_coverages=expand(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.antitargetcoverage.cnn", normal_id=get_normal_sample_ids())
    output:
        pooled_ref=f"{config['dirs']['pon_creation']}/pooled_reference.cnn"
    log:
        f"{config['dirs']['logs']}/cnvkit_reference_pooled/log.txt"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        """
        # Create pooled reference from clean normals, or fallback to flat reference
        if [ -s {input.normal_list} ]; then
            # Build reference from a panel of clean normals - need both target and antitarget
            echo "Building pooled reference from {input.normal_list}" > {log}
            
            # Extract clean normal IDs and create file patterns for both target and antitarget
            CLEAN_PATTERN=$(cat {input.normal_list} | sed 's/\\.targetcoverage\\.cnn//' | sed 's|^|{config[dirs][pon_creation]}/coverage/|' | sed 's|$|.*coverage.cnn|' | tr '\\n' ' ')
            
            echo "Using coverage files matching pattern: $CLEAN_PATTERN" &>> {log}
            cnvkit.py reference $CLEAN_PATTERN \\
                -f {config[reference_genome]} \\
                -o {output.pooled_ref} &>> {log}
        else
            # Fallback to a flat reference if no clean normals are available
            echo 'WARNING: No clean normals found. Creating a flat reference instead.' > {log}
            cnvkit.py reference \\
                -f {config[reference_genome]} \\
                -t {config[targets_bed]} \\
                -a {config[access_bed]} \\
                -o {output.pooled_ref} &>> {log}
        fi
        """
