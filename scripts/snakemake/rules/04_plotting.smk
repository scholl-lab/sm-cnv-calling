# ----------------------------------------------------------------------------------- #
# MODULE 4: PLOTTING AND VISUALIZATION
# ----------------------------------------------------------------------------------- #

rule cnvkit_scatter_genome:
    input:
        cnr=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cnr",
        cns=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cns",
        vcf=get_vcf_for_call,
    output:
        pdf=f"{config['dirs']['plots']}/scatter_genome/{{sample_id}}.genome.pdf"
    log:
        f"{config['dirs']['logs']}/cnvkit_scatter_genome/{{sample_id}}.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    resources:
        mem_mb=8000
    shell:
        """
        VCF_PARAM=""
        if [ -n "{input.vcf}" ]; then VCF_PARAM="-v {input.vcf}"; fi
        cnvkit.py scatter {input.cnr} -s {input.cns} $VCF_PARAM -o {output.pdf} &> {log}
        """

rule cnvkit_scatter_gene:
    input:
        cnr=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cnr",
        cns=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cns",
        vcf=get_vcf_for_call,
    output:
        pdf=f"{config['dirs']['plots']}/scatter_gene/{{sample_id}}.{{gene}}.pdf"
    log:
        f"{config['dirs']['logs']}/cnvkit_scatter_gene/{{sample_id}}.{{gene}}.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        """
        VCF_PARAM=""
        if [ -n "{input.vcf}" ]; then VCF_PARAM="-v {input.vcf}"; fi
        cnvkit.py scatter {input.cnr} -s {input.cns} $VCF_PARAM -g {wildcards.gene} -o {output.pdf} &> {log}
        """

rule cnvkit_diagram:
    input:
        cnr=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cnr",
        cns=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cns",
    output:
        pdf=f"{config['dirs']['plots']}/diagram/{{sample_id}}.diagram.pdf"
    log:
        f"{config['dirs']['logs']}/cnvkit_diagram/{{sample_id}}.log"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        "cnvkit.py diagram {input.cnr} -s {input.cns} -o {output.pdf} &> {log}"

rule cnvkit_heatmap_cohort:
    input:
        expand(f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cns", sample_id=get_tumor_samples())
    output:
        pdf=f"{config['dirs']['plots']}/heatmap_cohort.pdf"
    log:
        f"{config['dirs']['logs']}/cnvkit_heatmap_cohort/log.txt"
    conda:
        f"{config['conda_env_dir']}/cnvkit.yaml"
    shell:
        "cnvkit.py heatmap {input} -d -o {output.pdf} &> {log}"
