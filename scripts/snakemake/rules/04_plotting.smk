# ----------------------------------------------------------------------------------- #
# MODULE 4: PLOTTING AND VISUALIZATION
# ----------------------------------------------------------------------------------- #

rule cnvkit_scatter_genome:
    input:
        cnr=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cnr",
        cns=f"{config['dirs']['final_calls']}/{{sample_id}}.call.cns",
        vcf=get_vcf_for_call,
    output:
        pdf=f"{config['dirs']['plots']}/scatter_genome/{{sample_id}}.genome.pdf"
    params:
        tumor_sample_id_vcf=lambda w: SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"] if pd.notna(SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"]) and SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"] != "" else "",
        normal_sample_id_vcf=lambda w: SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"] if pd.notna(SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"]) and SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"] != "" else ""
    log:
        f"{config['dirs']['logs']}/cnvkit_scatter_genome/{{sample_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    resources:
        mem_mb=8000
    shell:
        """
        VCF_PARAM=""
        if [ -n "{input.vcf}" ]; then
            # Use tumor sample ID from VCF if specified, otherwise fallback to sample_id
            TUMOR_ID="{params.tumor_sample_id_vcf}"
            if [ -z "$TUMOR_ID" ]; then
                TUMOR_ID="{wildcards.sample_id}"
            fi
            VCF_PARAM="-v {input.vcf} -i $TUMOR_ID"
            echo "Using tumor sample ID for plotting: $TUMOR_ID" >> {log}
            
            # Add normal sample ID if provided
            if [ -n "{params.normal_sample_id_vcf}" ]; then
                VCF_PARAM="$VCF_PARAM -n {params.normal_sample_id_vcf}"
                echo "Using normal sample ID for plotting: {params.normal_sample_id_vcf}" >> {log}
            fi
        fi
        cnvkit.py scatter {input.cnr} -s {input.cns} $VCF_PARAM -o {output.pdf} &> {log}
        """

rule cnvkit_scatter_gene:
    input:
        cnr=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cnr",
        cns=f"{config['dirs']['final_calls']}/{{sample_id}}.call.cns",
        vcf=get_vcf_for_call,
    output:
        pdf=f"{config['dirs']['plots']}/scatter_gene/{{sample_id}}/{{gene}}.pdf"
    params:
        tumor_sample_id_vcf=lambda w: SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"] if pd.notna(SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"]) and SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"] != "" else "",
        normal_sample_id_vcf=lambda w: SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"] if pd.notna(SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"]) and SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"] != "" else ""
    log:
        f"{config['dirs']['logs']}/cnvkit_scatter_gene/{{sample_id}}.{{gene}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        """
        # Create sample-specific directory for gene plots
        mkdir -p $(dirname {output.pdf})
        
        VCF_PARAM=""
        if [ -n "{input.vcf}" ]; then
            # Use tumor sample ID from VCF if specified, otherwise fallback to sample_id
            TUMOR_ID="{params.tumor_sample_id_vcf}"
            if [ -z "$TUMOR_ID" ]; then
                TUMOR_ID="{wildcards.sample_id}"
            fi
            VCF_PARAM="-v {input.vcf} -i $TUMOR_ID"
            echo "Using tumor sample ID for gene plotting: $TUMOR_ID" >> {log}
            
            # Add normal sample ID if provided
            if [ -n "{params.normal_sample_id_vcf}" ]; then
                VCF_PARAM="$VCF_PARAM -n {params.normal_sample_id_vcf}"
                echo "Using normal sample ID for gene plotting: {params.normal_sample_id_vcf}" >> {log}
            fi
        fi
        cnvkit.py scatter {input.cnr} -s {input.cns} $VCF_PARAM -g {wildcards.gene} -o {output.pdf} &> {log}
        """

rule cnvkit_diagram:
    input:
        cnr=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cnr",
        cns=f"{config['dirs']['final_calls']}/{{sample_id}}.call.cns",
    output:
        pdf=f"{config['dirs']['plots']}/diagram/{{sample_id}}.diagram.pdf"
    log:
        f"{config['dirs']['logs']}/cnvkit_diagram/{{sample_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        "cnvkit.py diagram {input.cnr} -s {input.cns} -o {output.pdf} &> {log}"

rule cnvkit_heatmap_cohort:
    input:
        expand(f"{config['dirs']['final_calls']}/{{sample_id}}.call.cns", sample_id=get_tumor_samples())
    output:
        pdf=f"{config['dirs']['plots']}/heatmap_cohort.pdf"
    log:
        f"{config['dirs']['logs']}/cnvkit_heatmap_cohort/log.txt"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        "cnvkit.py heatmap {input} -d -o {output.pdf} &> {log}"
