# ----------------------------------------------------------------------------------- #
# MODULE 3: CNV CALLING AND EXPORT
# ----------------------------------------------------------------------------------- #

rule cnvkit_batch_tumor:
    input:
        tumor_bam=lambda w: SAMPLES.loc[w.sample_id, "tumor_bam"],
        pooled_ref=f"{config['dirs']['pon_creation']}/pooled_reference.cnn"
    output:
        cnr=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cnr",
        cns=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cns"
    log:
        f"{config['dirs']['logs']}/cnvkit_batch_tumor/{{sample_id}}.log"
    params:
        output_dir=config['dirs']['cnvkit_runs']
    conda:
        config["conda_envs"]["cnvkit"]
    threads: config["default_threads"]
    resources:
        mem_mb=config["cnvkit_batch_mem_mb"]
    shell:
        "cnvkit.py batch {input.tumor_bam} -r {input.pooled_ref} -p {threads} --exclude {config[blacklist]} --segment-threshold {config[segment_threshold]} -d {params.output_dir} &> {log}"

rule cnvkit_call:
    input:
        cns=f"{config['dirs']['cnvkit_runs']}/{{sample_id}}.cns",
        purity_file=f"{config['dirs']['purity_values']}/{{sample_id}}.purity.txt",
        vcf=get_vcf_for_call
    output:
        call_cns=f"{config['dirs']['final_calls']}/{{sample_id}}.call.cns"
    log:
        f"{config['dirs']['logs']}/cnvkit_call/{{sample_id}}.log"
    conda:
        config["conda_envs"]["cnvkit"]
    shell:
        """
        VCF_PARAM=""
        if [ -n "{input.vcf}" ]; then
            VCF_PARAM="-v {input.vcf}"
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
    shell:
        """
        cnvkit.py export vcf {input.call_cns} -i {params.sample_id} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
