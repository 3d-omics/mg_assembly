rule report__step__reads:
    """Collect all reports for the reads step"""
    input:
        [
            READS / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],
    output:
        REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    retries: 5
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
    """Collect all reports for the preprocess step"""
    input:
        fastp=[
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        fastqc=[
            FASTP / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],
        kraken2=[
            KRAKEN2 / kraken2_db / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken2_db in KRAKEN2_DBS
        ],
        bowtie2=[
            PRE_BOWTIE2 / host_name / f"{sample_id}.{library_id}.{report}"
            for host_name in HOST_NAMES
            for report in BAM_REPORTS
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        REPORT_STEP / "preprocess.html",
    log:
        REPORT_STEP / "preprocess.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    retries: 5
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


rule report__step__prokaryotes:
    """Collect all reports for the bowtie2 step when mapping to a mag catalogue"""
    input:
        reports=[
            QUANT_BOWTIE2 / f"{sample_id}.{library_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in ["stats.txt", "flagstats.txt"]
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
        rules.report__step__reads.output,
        rules.report__step__preprocess.output,
        rules.report__step__prokaryotes.output,
