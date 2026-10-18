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


if CNVKIT_AUTOBIN:

    rule cnvkit_autobin:
        """Calibrate target and antitarget bin sizes automatically from cohort BAMs."""
        input:
            bams=lambda wc: [
                BAM_PATH_MAP[bn] for bn in (get_normal_bams() or get_unique_bams()[:1])
            ],
            fasta=REF,
            targets=get_targets_bed_for_binning,
            access=get_access_bed,
            ref_flat=[REF_FLAT] if REF_FLAT else [],
        output:
            targets=os.path.join(CNVKIT_REF_DIR, "targets.bed"),
            antitargets=os.path.join(CNVKIT_REF_DIR, "antitargets.bed"),
        params:
            annotate_arg=f"--annotate {REF_FLAT}" if REF_FLAT else "",
            extra=CNVKIT_AUTOBIN_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wc, attempt: 8000 * attempt,
            runtime=lambda wc, attempt: min(240 * attempt, 480),
        conda:
            "../envs/cnvkit.yaml"
        log:
            os.path.join(LOG_DIR, "cnvkit_autobin/autobin.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/cnvkit_autobin/autobin.tsv")
        shell:
            r"""
            cnvkit.py autobin {input.bams} \
                -f {input.fasta} \
                -t {input.targets} \
                -g {input.access} \
                --short-names \
                {params.annotate_arg} \
                {params.extra} \
                --target-output-bed {output.targets} \
                --antitarget-output-bed {output.antitargets} &> {log}
            """

else:

    rule cnvkit_target:
        """Prepare split, short-named target bins with optional refFlat gene annotation."""
        input:
            targets=get_targets_bed_for_binning,
            ref_flat=[REF_FLAT] if REF_FLAT else [],
        output:
            targets=os.path.join(CNVKIT_REF_DIR, "targets.bed"),
        params:
            annotate_arg=f"--annotate {REF_FLAT}" if REF_FLAT else "",
            avg_size_arg=(
                f"--avg-size {CNVKIT_TARGET_AVG_SIZE}" if CNVKIT_TARGET_AVG_SIZE > 0 else ""
            ),
            extra=CNVKIT_TARGET_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wc, attempt: 4000 * attempt,
            runtime=lambda wc, attempt: min(60 * attempt, 480),
        conda:
            "../envs/cnvkit.yaml"
        log:
            os.path.join(LOG_DIR, "cnvkit_target/target.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/cnvkit_target/target.tsv")
        shell:
            r"""
            cnvkit.py target {input.targets} \
                --split \
                --short-names \
                {params.annotate_arg} \
                {params.avg_size_arg} \
                {params.extra} \
                -o {output.targets} &> {log}
            """

    rule cnvkit_antitarget:
        """Generate off-target (antitarget) bins from accessible genome regions."""
        input:
            targets=os.path.join(CNVKIT_REF_DIR, "targets.bed"),
            access=get_access_bed,
        output:
            antitargets=os.path.join(CNVKIT_REF_DIR, "antitargets.bed"),
        params:
            avg_size_arg=(
                f"--avg-size {CNVKIT_ANTITARGET_AVG_SIZE}"
                if CNVKIT_ANTITARGET_AVG_SIZE > 0
                else ""
            ),
            min_size_arg=(
                f"--min-size {CNVKIT_ANTITARGET_MIN_SIZE}"
                if CNVKIT_ANTITARGET_MIN_SIZE > 0
                else ""
            ),
            extra=CNVKIT_ANTITARGET_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wc, attempt: 4000 * attempt,
            runtime=lambda wc, attempt: min(60 * attempt, 480),
        conda:
            "../envs/cnvkit.yaml"
        log:
            os.path.join(LOG_DIR, "cnvkit_antitarget/antitarget.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/cnvkit_antitarget/antitarget.tsv")
        shell:
            r"""
            cnvkit.py antitarget {input.targets} \
                -g {input.access} \
                {params.avg_size_arg} \
                {params.min_size_arg} \
                {params.extra} \
                -o {output.antitargets} &> {log}
            """


rule cnvkit_coverage_target:
    """Compute on-target read-depth coverage once per unique BAM."""
    input:
        bam=get_bam_by_basename,
        targets=os.path.join(CNVKIT_REF_DIR, "targets.bed"),
    output:
        target=os.path.join(CNVKIT_COV_DIR, "{bam_basename}.targetcoverage.cnn"),
    params:
        extra=CNVKIT_COVERAGE_EXTRA,
    threads: 4
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 8000 * attempt,
        runtime=lambda wc, attempt: min(240 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_coverage_target/{bam_basename}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_coverage_target/{bam_basename}.tsv")
    shell:
        r"""
        cnvkit.py coverage {input.bam} {input.targets} \
            -p {threads} \
            {params.extra} \
            -o {output.target} &> {log}
        """


rule cnvkit_coverage_antitarget:
    """Compute off-target read-depth coverage once per unique BAM."""
    input:
        bam=get_bam_by_basename,
        antitargets=os.path.join(CNVKIT_REF_DIR, "antitargets.bed"),
    output:
        antitarget=os.path.join(CNVKIT_COV_DIR, "{bam_basename}.antitargetcoverage.cnn"),
    params:
        extra=CNVKIT_COVERAGE_EXTRA,
    threads: 4
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 8000 * attempt,
        runtime=lambda wc, attempt: min(240 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_coverage_antitarget/{bam_basename}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_coverage_antitarget/{bam_basename}.tsv")
    shell:
        r"""
        cnvkit.py coverage {input.bam} {input.antitargets} \
            -p {threads} \
            {params.extra} \
            -o {output.antitarget} &> {log}
        """


rule cnvkit_fix:
    """Normalize tumor target and antitarget coverage against the PoN reference."""
    input:
        target=lambda wc: os.path.join(
            CNVKIT_COV_DIR,
            f"{get_unit_tumor_bam_basename(wc.unit)}.targetcoverage.cnn",
        ),
        antitarget=lambda wc: os.path.join(
            CNVKIT_COV_DIR,
            f"{get_unit_tumor_bam_basename(wc.unit)}.antitargetcoverage.cnn",
        ),
        reference=get_cnvkit_reference,
    output:
        cnr=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cnr"),
    params:
        extra=CNVKIT_FIX_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(120 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_fix/{unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_fix/{unit}.tsv")
    shell:
        r"""
        cnvkit.py fix {input.target} {input.antitarget} {input.reference} \
            -i {wildcards.unit} \
            {params.extra} \
            -o {output.cnr} &> {log}
        """


rule cnvkit_segment:
    """Segment normalized log2 copy ratios with optional allele-aware VCF BAF."""
    input:
        cnr=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cnr"),
        vcf=get_unit_vcf_input,
    output:
        cns=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cns"),
    params:
        method=CNVKIT_SEGMENT_METHOD,
        threshold=CNVKIT_SEGMENT_THRESHOLD,
        drop_low_coverage_arg="--drop-low-coverage" if CNVKIT_DROP_LOW_COVERAGE else "",
        vcf_args=get_unit_vcf_args,
        extra=CNVKIT_SEGMENT_EXTRA,
    threads: 4
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 8000 * attempt,
        runtime=lambda wc, attempt: min(240 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_segment/{unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_segment/{unit}.tsv")
    shell:
        r"""
        cnvkit.py segment {input.cnr} \
            -m {params.method} \
            -t {params.threshold} \
            {params.drop_low_coverage_arg} \
            {params.vcf_args} \
            -p {threads} \
            {params.extra} \
            -o {output.cns} &> {log}
        """


rule cnvkit_bintest:
    """Call focal bin-level amplifications and deletions."""
    input:
        cnr=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cnr"),
        cns=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cns"),
    output:
        bintest_cns=os.path.join(CNVKIT_CALLS_DIR, "{unit}.bintest.cns"),
    params:
        extra=CNVKIT_BINTEST_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(120 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_bintest/{unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_bintest/{unit}.tsv")
    shell:
        r"""
        cnvkit.py bintest {input.cnr} \
            -s {input.cns} \
            {params.extra} \
            -o {output.bintest_cns} &> {log}
        """


