# ----------------------------------------------------------------------------------- #
# MODULE 2: PANEL OF NORMALS (PoN) GENERATION
# ----------------------------------------------------------------------------------- #

ruleorder: cnvkit_reference_pooled > cnvkit_batch_tumor

# Prepare annotated targets BED file (needed before antitarget creation)
rule cnvkit_target:
    input:
        targets_raw=config["targets_bed"]
    output:
        targets_annotated=f"{config['dirs']['pon_creation']}/targets_annotated.bed"
    params:
        annotate_param=f"--annotate {config['annotate_refFlat']}" if config.get("annotate_refFlat") else ""
    log:
        f"{config['dirs']['logs']}/cnvkit_target/log.txt"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        """
        if [ -n "{params.annotate_param}" ]; then
            echo "Annotating targets with gene information from {config[annotate_refFlat]}" > {log}
            cnvkit.py target {input.targets_raw} {params.annotate_param} --split -o {output.targets_annotated} &>> {log}
        else
            echo "No gene annotation file provided, using targets as-is" > {log}
            cp {input.targets_raw} {output.targets_annotated} &>> {log}
        fi
        """

# Global rule to create antitarget regions (needed by both PoN and PureCN workflows)
rule create_antitarget_bed:
    input:
        targets=rules.cnvkit_target.output.targets_annotated,
        access=config["access_bed"]
    output:
        antitargets=f"{config['dirs']['pon_creation']}/antitargets.bed"
    log:
        f"{config['dirs']['logs']}/create_antitarget_bed/log.txt"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        "cnvkit.py antitarget {input.targets} -g {input.access} -o {output.antitargets} &> {log}"

rule cnvkit_coverage_normals_target:
    input:
        bam=lambda w: get_normal_bams_for_pon()[w.normal_id],
        targets=rules.cnvkit_target.output.targets_annotated
    output:
        target=f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.targetcoverage.cnn"
    log:
        f"{config['dirs']['logs']}/cnvkit_coverage_normals_target/{{normal_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    threads: config["default_threads"]
    shell:
        "cnvkit.py coverage {input.bam} {input.targets} -p {threads} -o {output.target} &> {log}"

rule cnvkit_coverage_normals_antitarget:
    input:
        bam=lambda w: get_normal_bams_for_pon()[w.normal_id],
        antitargets=rules.create_antitarget_bed.output.antitargets
    output:
        antitarget=f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.antitargetcoverage.cnn"
    log:
        f"{config['dirs']['logs']}/cnvkit_coverage_normals_antitarget/{{normal_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    threads: config["default_threads"]
    shell:
        "cnvkit.py coverage {input.bam} {input.antitargets} -p {threads} -o {output.antitarget} &> {log}"

rule cnvkit_metrics_normals:
    input:
        target_coverages=expand(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.targetcoverage.cnn", normal_id=get_normal_sample_ids()),
        antitarget_coverages=expand(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.antitargetcoverage.cnn", normal_id=get_normal_sample_ids())
    output:
        metrics_file=f"{config['dirs']['pon_creation']}/qc/all_normals.metrics.txt"
    log:
        f"{config['dirs']['logs']}/cnvkit_metrics_normals/log.txt"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        "cnvkit.py metrics {input.target_coverages} > {output.metrics_file} 2> {log}"

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
        antitarget_coverages=expand(f"{config['dirs']['pon_creation']}/coverage/{{normal_id}}.antitargetcoverage.cnn", normal_id=get_normal_sample_ids()),
        targets_annotated=rules.cnvkit_target.output.targets_annotated,
        antitargets=rules.create_antitarget_bed.output.antitargets
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
            
            # Extract clean normal IDs and create explicit file lists for both target and antitarget
            CLEAN_TARGETS=$(cat {input.normal_list} | tr '\n' ' ')
            CLEAN_ANTITARGETS=$(cat {input.normal_list} | sed 's/\.targetcoverage\.cnn/.antitargetcoverage.cnn/' | tr '\n' ' ')
            CLEAN_FILES="$CLEAN_TARGETS $CLEAN_ANTITARGETS"
            
            echo "Using clean normal files: $CLEAN_FILES" &>> {log}
            cnvkit.py reference $CLEAN_FILES \\
                -f {config[reference_genome]} \\
                -o {output.pooled_ref} &>> {log}
        else
            # Fallback to a flat reference if no clean normals are available
            echo 'WARNING: No clean normals found. Creating a flat reference instead.' > {log}
            cnvkit.py reference \\
                -f {config[reference_genome]} \\
                -t {input.targets_annotated} \\
                -a {input.antitargets} \\
                -o {output.pooled_ref} &>> {log}
        fi
        """
