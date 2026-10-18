"""CNVkit rules: access, target/antitarget binning, PoN reference, fix, segment, call, exports."""


rule cnvkit_access:
    """Compute or filter mappable accessible regions BED for GRCh38."""
    input:
        fasta=REF,
        access_raw=[REF_ACCESS_BED] if REF_ACCESS_BED else [],
        blacklist=[REF_BLACKLIST_BED] if REF_BLACKLIST_BED else [],
    output:
        access=os.path.join(CNVKIT_REF_DIR, "access.bed"),
    params:
        min_gap_size=CNVKIT_ACCESS_MIN_GAP_SIZE,
        blacklist_arg=f"-x {REF_BLACKLIST_BED}" if REF_BLACKLIST_BED else "",
        extra=CNVKIT_ACCESS_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(120 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_access/access.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_access/access.tsv")
    shell:
        r"""
        if [ -n "{input.access_raw}" ] && [ -n "{input.blacklist}" ]; then
            bedtools subtract -a {input.access_raw} -b {input.blacklist} > {output.access} 2> {log}
        elif [ -n "{input.access_raw}" ]; then
            cp {input.access_raw} {output.access} 2> {log}
        else
            cnvkit.py access {input.fasta} \
                -s {params.min_gap_size} \
                {params.blacklist_arg} \
                {params.extra} \
                -o {output.access} &> {log}
        fi
        """


rule cnvkit_filter_targets:
    """Subtract genomic blacklist intervals from capture target BED before binning."""
    input:
        targets=TARGETS_BED,
        blacklist=[REF_BLACKLIST_BED] if REF_BLACKLIST_BED else [],
    output:
        targets_clean=os.path.join(CNVKIT_REF_DIR, "targets.clean.bed"),
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(60 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_filter_targets/filter_targets.log"),
    shell:
        r"""
        if [ -n "{input.blacklist}" ]; then
            bedtools subtract -a {input.targets} -b {input.blacklist} > {output.targets_clean} 2> {log}
        else
            cp {input.targets} {output.targets_clean} 2> {log}
        fi
        """


