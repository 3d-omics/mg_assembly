rule report_step_reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads_eval_fastqc.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_STEP,
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


rule report_step_pre:
    """Collect all reports for the preprocessing step"""
    input:
        rules.pre_eval_fastp.input,
        rules.pre_eval_fastp_fastqc.input,
        rules.pre_eval_samtools.input,
        rules.pre_eval_nonhost_fastqc.input,
        rules.pre_eval_kraken2.input,
    output:
        html=REPORT_STEP / "preprocessing.html",
    log:
        REPORT_STEP / "preprocessing.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title preprocessing \
            --force \
            --filename preprocessing \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report_step_assemble:
    """Collect all reports for the assemble step"""
    input:
        rules.assemble_eval_quast.input,
        rules.assemble_eval_samtools.input,
    output:
        html=REPORT_STEP / "assemble.html",
    log:
        REPORT_STEP / "assemble.log",
    params:
        dir=REPORT_STEP,
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


rule report_step:
    input:
        REPORT_STEP / "reads.html",
        REPORT_STEP / "preprocessing.html",
        REPORT_STEP / "assemble.html",
