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


rule cnvkit_sex:
    """Infer chromosomal sex from target and antitarget coverage per BAM."""
    input:
        target=os.path.join(CNVKIT_COV_DIR, "{bam_basename}.targetcoverage.cnn"),
        antitarget=os.path.join(CNVKIT_COV_DIR, "{bam_basename}.antitargetcoverage.cnn"),
    output:
        sex_tsv=os.path.join(CNVKIT_COV_DIR, "{bam_basename}.sex.tsv"),
    params:
        male_ref_arg="-y" if CNVKIT_MALE_REFERENCE else "",
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(60 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_sex/{bam_basename}.log"),
    shell:
        r"""
        cnvkit.py sex {input.target} {input.antitarget} \
            {params.male_ref_arg} \
            -o {output.sex_tsv} &> {log}
        """


rule resolve_unit_sex:
    """Resolve chromosomal sex for a CNV unit (sample sheet override or cnvkit sex)."""
    localrule: True
    input:
        sex_tsv=lambda wc: os.path.join(
            CNVKIT_COV_DIR, f"{get_unit_tumor_bam_basename(wc.unit)}.sex.tsv"
        ),
    output:
        sex_txt=os.path.join(CNVKIT_CALLS_DIR, "{unit}.sex.txt"),
    params:
        sex_override=lambda wc: CNV_UNITS[wc.unit]["sex"],
        rules_dir=RULES_DIR,
    log:
        os.path.join(LOG_DIR, "resolve_unit_sex/{unit}.log"),
    shell:
        r"""
        python3 -c '
import sys
from pathlib import Path
sys.path.insert(0, "{params.rules_dir}")
from helpers import parse_cnvkit_sex
sex = parse_cnvkit_sex("{input.sex_tsv}", "{params.sex_override}")
Path("{output.sex_txt}").parent.mkdir(parents=True, exist_ok=True)
Path("{output.sex_txt}").write_text(sex + "\n", encoding="utf-8")
' &> {log}
        """


rule cnvkit_metrics_normals:
    """Compute noise metrics across candidate PoN normal target coverages."""
    input:
        target_coverages=lambda wc: expand(
            os.path.join(CNVKIT_COV_DIR, "{bn}.targetcoverage.cnn"),
            bn=get_normal_bams(),
        ),
        antitarget_coverages=lambda wc: expand(
            os.path.join(CNVKIT_COV_DIR, "{bn}.antitargetcoverage.cnn"),
            bn=get_normal_bams(),
        ),
    output:
        metrics_file=os.path.join(CNVKIT_PON_QC_DIR, "all_normals.metrics.txt"),
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(120 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_metrics_normals/metrics.log"),
    shell:
        r"""
        if [ -n "{input.target_coverages}" ]; then
            cnvkit.py metrics {input.target_coverages} -o {output.metrics_file} &> {log}
        else
            printf "sample\tsegments\tstdev\tmad\tiqr\tbivar\n" > {output.metrics_file}
            echo "No normal BAMs in cohort; wrote empty normal metrics header." > {log}
        fi
        """


rule identify_clean_normals:
    """Filter normal coverages by bivar threshold and write verified .cnn paths."""
    localrule: True
    input:
        metrics=os.path.join(CNVKIT_PON_QC_DIR, "all_normals.metrics.txt"),
        target_coverages=lambda wc: expand(
            os.path.join(CNVKIT_COV_DIR, "{bn}.targetcoverage.cnn"),
            bn=get_normal_bams(),
        ),
        antitarget_coverages=lambda wc: expand(
            os.path.join(CNVKIT_COV_DIR, "{bn}.antitargetcoverage.cnn"),
            bn=get_normal_bams(),
        ),
    output:
        clean_list=os.path.join(CNVKIT_PON_QC_DIR, "clean_normals.list"),
    params:
        threshold=CNVKIT_PON_BIVAR_THRESHOLD,
        cov_dir=CNVKIT_COV_DIR,
        rules_dir=RULES_DIR,
    log:
        os.path.join(LOG_DIR, "identify_clean_normals/clean_normals.log"),
    shell:
        r"""
        python3 -c '
import sys
sys.path.insert(0, "{params.rules_dir}")
from helpers import identify_clean_normals
identify_clean_normals("{input.metrics}", "{params.cov_dir}", "{output.clean_list}", float({params.threshold}))
' &> {log}
        """


rule cnvkit_reference_pooled:
    """Build pooled normal reference from clean normals or fallback to flat reference."""
    input:
        normal_list=lambda wc: (
            [os.path.join(CNVKIT_PON_QC_DIR, "clean_normals.list")] if get_normal_bams() else []
        ),
        target_coverages=lambda wc: expand(
            os.path.join(CNVKIT_COV_DIR, "{bn}.targetcoverage.cnn"),
            bn=get_normal_bams(),
        ),
        antitarget_coverages=lambda wc: expand(
            os.path.join(CNVKIT_COV_DIR, "{bn}.antitargetcoverage.cnn"),
            bn=get_normal_bams(),
        ),
        fasta=REF,
        targets=os.path.join(CNVKIT_REF_DIR, "targets.bed"),
        antitargets=os.path.join(CNVKIT_REF_DIR, "antitargets.bed"),
    output:
        pooled_ref=os.path.join(CNVKIT_REF_DIR, "pooled_reference.cnn"),
    params:
        male_ref_arg="-y" if CNVKIT_MALE_REFERENCE else "",
        extra=CNVKIT_REFERENCE_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 8000 * attempt,
        runtime=lambda wc, attempt: min(240 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_reference_pooled/reference.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_reference_pooled/reference.tsv")
    shell:
        r"""
        if [ -n "{input.normal_list}" ] && [ -s "{input.normal_list}" ]; then
            CLEAN_FILES=$(tr '\n' ' ' < {input.normal_list})
            cnvkit.py reference $CLEAN_FILES \
                -f {input.fasta} \
                {params.male_ref_arg} \
                {params.extra} \
                -o {output.pooled_ref} &> {log}
        else
            echo "WARNING: No clean normals available; building flat reference." > {log}
            cnvkit.py reference \
                -f {input.fasta} \
                -t {input.targets} \
                -a {input.antitargets} \
                {params.male_ref_arg} \
                {params.extra} \
                -o {output.pooled_ref} &>> {log}
        fi
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


rule cnvkit_segmetrics:
    """Compute segment-level confidence intervals and dispersion statistics."""
    input:
        cnr=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cnr"),
        cns=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cns"),
    output:
        segmetrics_cns=os.path.join(CNVKIT_RUNS_DIR, "{unit}.segmetrics.cns"),
    params:
        extra=CNVKIT_SEGMETRICS_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(120 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_segmetrics/{unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_segmetrics/{unit}.tsv")
    shell:
        r"""
        cnvkit.py segmetrics {input.cnr} \
            -s {input.cns} \
            --ci --pi --sem --iqr \
            {params.extra} \
            -o {output.segmetrics_cns} &> {log}
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


rule cnvkit_call:
    """Call integer copy numbers and allele-specific CN using purity, ploidy, sex, and BAF."""
    input:
        cns=os.path.join(CNVKIT_RUNS_DIR, "{unit}.segmetrics.cns"),
        purity_ploidy=os.path.join(PURITY_DIR, "{unit}.purity_ploidy.tsv"),
        sex_txt=os.path.join(CNVKIT_CALLS_DIR, "{unit}.sex.txt"),
        vcf=get_unit_vcf_input,
    output:
        call_cns=os.path.join(CNVKIT_CALLS_DIR, "{unit}.call.cns"),
    params:
        male_reference=CNVKIT_MALE_REFERENCE,
        call_method=CNVKIT_CALL_METHOD,
        call_filter=CNVKIT_CALL_FILTER,
        drop_low_coverage=CNVKIT_DROP_LOW_COVERAGE,
        vcf_args=get_unit_vcf_args,
        extra=CNVKIT_CALL_EXTRA,
        rules_dir=RULES_DIR,
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 8000 * attempt,
        runtime=lambda wc, attempt: min(120 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_call/{unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/cnvkit_call/{unit}.tsv")
    shell:
        r"""
        CALL_ARGS=$(python3 -c '
import sys
from pathlib import Path
sys.path.insert(0, "{params.rules_dir}")
from helpers import build_cnvkit_call_args
sex = Path("{input.sex_txt}").read_text(encoding="utf-8").strip()
print(build_cnvkit_call_args(
    purity_ploidy_tsv="{input.purity_ploidy}",
    sex=sex,
    male_reference=bool({params.male_reference}),
    call_method="{params.call_method}",
    call_filter="{params.call_filter}",
    drop_low_coverage=bool({params.drop_low_coverage}),
    vcf_args={params.vcf_args!r},
    extra={params.extra!r},
))
')
        echo "Running: cnvkit.py call {input.cns} $CALL_ARGS -o {output.call_cns}" > {log}
        cnvkit.py call {input.cns} $CALL_ARGS -o {output.call_cns} &>> {log}
        """


rule cnvkit_export_seg:
    """Export segmentation in DNAcopy SEG format for PureCN, IGV, and GISTIC2."""
    input:
        cns=os.path.join(CNVKIT_RUNS_DIR, "{unit}.cns"),
    output:
        seg=os.path.join(CNVKIT_EXPORTS_DIR, "{unit}.seg"),
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(60 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_export_seg/{unit}.log"),
    shell:
        r"""
        cnvkit.py export seg {input.cns} -o {output.seg} &> {log}
        """


rule cnvkit_export_bed:
    """Export non-neutral integer copy-number segments in BED format."""
    input:
        call_cns=os.path.join(CNVKIT_CALLS_DIR, "{unit}.call.cns"),
        sex_txt=os.path.join(CNVKIT_CALLS_DIR, "{unit}.sex.txt"),
    output:
        bed=os.path.join(CNVKIT_EXPORTS_DIR, "{unit}.cnv.bed"),
    params:
        male_ref_arg="-y" if CNVKIT_MALE_REFERENCE else "",
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(60 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_export_bed/{unit}.log"),
    shell:
        r"""
        SAMPLE_SEX=$(cat {input.sex_txt})
        cnvkit.py export bed {input.call_cns} \
            --show variant \
            -i {wildcards.unit} \
            {params.male_ref_arg} \
            -x "$SAMPLE_SEX" \
            -o {output.bed} &> {log}
        """


rule cnvkit_export_vcf:
    """Export final copy-number calls as bgzipped and tabix-indexed VCF."""
    input:
        call_cns=os.path.join(CNVKIT_CALLS_DIR, "{unit}.call.cns"),
        sex_txt=os.path.join(CNVKIT_CALLS_DIR, "{unit}.sex.txt"),
    output:
        vcf=os.path.join(CNVKIT_EXPORTS_DIR, "{unit}.cnv.vcf.gz"),
        tbi=os.path.join(CNVKIT_EXPORTS_DIR, "{unit}.cnv.vcf.gz.tbi"),
    params:
        male_ref_arg="-y" if CNVKIT_MALE_REFERENCE else "",
    retries: 2
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: min(60 * attempt, 480),
    conda:
        "../envs/cnvkit.yaml"
    log:
        os.path.join(LOG_DIR, "cnvkit_export_vcf/{unit}.log"),
    shell:
        r"""
        SAMPLE_SEX=$(cat {input.sex_txt})
        cnvkit.py export vcf {input.call_cns} \
            -i {wildcards.unit} \
            {params.male_ref_arg} \
            -x "$SAMPLE_SEX" \
            2> {log} | bgzip -c > {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf} &>> {log}
        """
