rule report__step__prokaryotes:
    """Collect all reports for the bowtie2 step when mapping to a mag catalogue"""
    input:
        reports=[
            QUANT_BOWTIE2
            / f"drep.{secondary_ani}"
            / f"{sample_id}.{library_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in ["stats.tsv", "flagstats.txt"]
            for secondary_ani in SECONDARY_ANIS
        ],
    output:
        REPORT_STEP / "prokaryotes.html",
    log:
        REPORT_STEP / "prokaryotes.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    retries: 5
    shell:
        """
        multiqc \
            --title prokaryotes \
            --force \
            --filename prokaryotes \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input.reports} \
        2> {log} 1>&2
        """


rule report__step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report__step__prokaryotes.output,
