rule report__sample__multiqc__:
    input:
        reads_fastqc=lambda w: [
            READS / f"{w.sample_id}.{w.library_id}_{end}_fastqc.zip"
            for end in [1, 2]
        ], 
        pre_fastp_fastqc=lambda w: [
            FASTP / f"{w.sample_id}.{w.library_id}_{end}_fastqc.zip"
            for end in [1, 2]
        ],
        pre_fastp_html=FASTP / "{sample_id}.{library_id}_fastp.json",
        bowtie2=lambda w: [
            PRE_BOWTIE2 / host_name / f"{w.sample_id}.{w.library_id}.{report}"
            for host_name in HOST_NAMES
            for report in BAM_REPORTS
        ],
        kraken2=lambda w: [
            KRAKEN2 / kraken2_db / f"{w.sample_id}.{w.library_id}.report"
            for kraken2_db in KRAKEN2_DBS
        ],
        quantify_bowtie2=lambda w: [
            PROK_QUANT / f"{w.sample_id}.{w.library_id}.{extension}"
            for extension in ["stats.txt", "flagstats.txt"]
        ],
    output:
        html=REPORT_SAMPLE / "{sample_id}.{library_id}.html",
    log:
        REPORT_SAMPLE / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_SAMPLE,
        filename=lambda wildcards: f"{wildcards.sample_id}.{wildcards.library_id}",
    shell:
        """
        multiqc \
            --filename {params.filename} \
            --title {params.filename} \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__sample:
    input:
        [
            REPORT_SAMPLE / f"{sample_id}.{library_id}.html"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
