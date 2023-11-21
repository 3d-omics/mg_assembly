rule report_step_reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads_fastqc.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        multiqc \
            --filename reads \
            --title reads \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__preprocess:
    """Collect all reports for the preprocessing step"""
    input:
        rules.preprocess__fastqc.input,
        rules.preprocess__fastp__json.input,
        rules.preprocess__samtools.input,
        rules.preprocess__kraken2.input,
    output:
        html=REPORT_STEP / "preprocess.html",
    log:
        REPORT_STEP / "preprocess.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        multiqc \
            --title preprocess \
            --force \
            --filename preprocess \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report__step__assemble:
    """Collect all reports for the assemble step"""
    input:
        rules.assemble__quast.input,
        rules.assemble__samtools.input,
    output:
        html=REPORT_STEP / "assemble.html",
    log:
        REPORT_STEP / "assemble.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        multiqc \
            --title assemble \
            --force \
            --filename assemble \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report__step__bin:
    """Collect all reports for the bin step"""
    input:
        rules.bin__quast.input,
    output:
        html=REPORT_STEP / "bin.html",
    log:
        REPORT_STEP / "bin.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        multiqc \
            --title bin \
            --force \
            --filename bin \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report__step__dereplicate:
    """Collect all reports for the dereplicate step"""
    input:
        rules.dereplicate__quast.input,
        rules.dereplicate__samtools.input,
    output:
        html=REPORT_STEP / "dereplicate.html",
    log:
        REPORT_STEP / "dereplicate.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        multiqc \
            --title dereplicate \
            --force \
            --filename dereplicate \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report_step:
    """Report for all steps"""
    input:
        REPORT_STEP / "reads.html",
        REPORT_STEP / "preprocess.html",
        REPORT_STEP / "assemble.html",
        REPORT_STEP / "bin.html",


rule report__step__with_dereplicate:
    """Report all steps + dereplicate"""
    input:
        rules.report_step.input,
        REPORT_STEP / "dereplicate.html",


localrules:
    report_step_reads,
    report__step__preprocess,
    report__step__assemble,
    report__step__bin,
    report_step,
    report__step__dereplicate,
    report__step__with_dereplicate,
