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
    conda:
        "report.yml"
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


rule report_step_dereplicate:
    """Collect all reports for the dereplicate step"""
    input:
        rules.dereplicate_eval_quast.input,
        rules.dereplicate_eval_samtools.input,
    output:
        html=REPORT_STEP / "dereplicate.html",
    log:
        REPORT_STEP / "dereplicate.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_STEP,
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
    input:
        REPORT_STEP / "reads.html",
        REPORT_STEP / "preprocessing.html",
        REPORT_STEP / "assemble.html",
        # REPORT_STEP / "bin.html",
        REPORT_STEP / "metabin.html",


rule report_step_with_dereplicate:
    input:
        rules.report_step.input,
        REPORT_STEP / "dereplicate.html",


localrules:
    report_step_reads,
    report_step_pre,
    report_step_assemble,
    report_step,
    report_step_with_dereplicate,
