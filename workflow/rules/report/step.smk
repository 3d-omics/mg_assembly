rule report__step__reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads__fastqc.input,
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
    """Collcect all reports from the assemble step"""
    input:
        QUAST,
    output:
        REPORT_STEP / "assemble.html",
    log:
        "_env.yml",
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
            {input} \
        2> {log} 1>&2
        """


rule report__step:
    """Report for all steps"""
    input:
        REPORT_STEP / "reads.html",
        REPORT_STEP / "preprocess.html",
        REPORT_STEP / "assemble.html",
        # REPORT_STEP / "quantify.html",
