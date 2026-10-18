# ----------------------------------------------------------------------------------- #
# MODULE 3: CNV CALLING AND EXPORT
# ----------------------------------------------------------------------------------- #

rule cnvkit_batch_tumor:
    input:
        tumor_bam=lambda w: SAMPLES.loc[w.sample_id, "tumor_bam"],
        pooled_ref=f"{config['dirs']['pon_creation']}/pooled_reference.cnn"
    output:
        cnr=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cnr",
        cns=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cns",
        targetcoverage=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.targetcoverage.cnn",
        antitargetcoverage=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.antitargetcoverage.cnn"
    log:
        f"{config['dirs']['logs']}/cnvkit_batch_tumor/{{sample_id}}.log"
    params:
        output_dir=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}_batch",
        bam_basename=lambda w: get_cnvkit_basename(w.sample_id)
    conda:
        config["conda_envs"]["cnvkit"]
    threads: config["default_threads"]
    resources:
        mem_mb=config["cnvkit_batch_mem_mb"]
    shell:
        """
        # Create sample-specific output directory to avoid conflicts
        mkdir -p {params.output_dir}
        
        # Run CNVkit batch - output files will be named based on BAM basename
        # Note: --annotate flag cannot be used with -r/--reference, annotations should be in the reference
        cnvkit.py batch {input.tumor_bam} -r {input.pooled_ref} -p {threads} -d {params.output_dir} &> {log}
        
        # Move all CNVkit output files to expected sample_id-based names in main output directory
        mv {params.output_dir}/{params.bam_basename}.cnr {output.cnr}
        mv {params.output_dir}/{params.bam_basename}.cns {output.cns}
        mv {params.output_dir}/{params.bam_basename}.targetcoverage.cnn {output.targetcoverage}
        mv {params.output_dir}/{params.bam_basename}.antitargetcoverage.cnn {output.antitargetcoverage}
        
        # Handle additional files that may be created (bintest.cns, call.cns) but are optional
        if [ -f "{params.output_dir}/{params.bam_basename}.bintest.cns" ]; then
            mv {params.output_dir}/{params.bam_basename}.bintest.cns {config[dirs][cnvkit_runs]}/{wildcards.sample_id}.bintest.cns
        fi
        if [ -f "{params.output_dir}/{params.bam_basename}.call.cns" ]; then
            mv {params.output_dir}/{params.bam_basename}.call.cns {config[dirs][cnvkit_runs]}/{wildcards.sample_id}.call.cns
        fi
        
        # Clean up sample-specific directory after moving files
        rm -rf {params.output_dir}
        """

rule cnvkit_call:
    input:
        cns=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cns",
        purity_file=f"{config['dirs']['purity_values']}/{{sample_id}}.purity.txt",
        vcf=get_vcf_for_call
    output:
        call_cns=f"{config['dirs']['final_calls']}/{{sample_id}}.call.cns"
    params:
        tumor_sample_id_vcf=lambda w: SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"] if pd.notna(SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"]) and SAMPLES.loc[w.sample_id, "tumor_sample_id_vcf"] != "" else "",
        normal_sample_id_vcf=lambda w: SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"] if pd.notna(SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"]) and SAMPLES.loc[w.sample_id, "normal_sample_id_vcf"] != "" else ""
    log:
        f"{config['dirs']['logs']}/cnvkit_call/{{sample_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    resources:
        mem_mb=config["cnvkit_call_mem_mb"]
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
            echo "Using tumor sample ID: $TUMOR_ID" >> {log}
            
            # Add normal sample ID if provided
            if [ -n "{params.normal_sample_id_vcf}" ]; then
                VCF_PARAM="$VCF_PARAM -n {params.normal_sample_id_vcf}"
                echo "Using normal sample ID: {params.normal_sample_id_vcf}" >> {log}
            fi
            
            echo "VCF parameters: $VCF_PARAM" >> {log}
        fi
        cnvkit.py call {input.cns} $VCF_PARAM --purity $(cat {input.purity_file}) -m clonal -o {output.call_cns} &> {log}
        """

rule cnvkit_export_vcf:
    input:
        call_cns=f"{config['dirs']['final_calls']}/{{sample_id}}.call.cns"
    output:
        vcf=f"{config['dirs']['final_vcfs']}/{{sample_id}}.cnv.vcf.gz",
        tbi=f"{config['dirs']['final_vcfs']}/{{sample_id}}.cnv.vcf.gz.tbi"
    params:
        sample_id="{sample_id}"
    log:
        f"{config['dirs']['logs']}/cnvkit_export_vcf/{{sample_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    resources:
        mem_mb=config["cnvkit_export_mem_mb"]
    shell:
        """
        cnvkit.py export vcf {input.call_cns} -i {params.sample_id} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
