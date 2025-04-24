rule preprocess__multiqc:
    input:
        reads=[
            PRE_READS / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
        fastp=[
            PRE_FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        bowtie2=[
            PRE_BOWTIE2 / f"{host}.{sample_id}.{library_id}.{report}"
            for host in HOST_NAMES
            for report in BAM_REPORTS
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        clean=[
            PRE_CLEAN / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],
        nonpareil=[
            PRE_NONPAREIL / f"{sample_id}.{library_id}.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        kraken2=[
            PRE_BRACKEN / kraken2_db / f"{sample_id}.k2report"
            for sample_id in SAMPLES
            for kraken2_db in KRAKEN2_DBS
        ],
    output:
        RESULTS / "preprocess.html",
        RESULTS / "preprocess_data.zip",
    log:
        RESULTS / "preprocess.log",
    params:
        extra="--title preprocess --dirs --fullnames --fn_as_s_name --force",
    resources:
        mem_mb=double_ram(4 * 1024),
        runtime=6 * 60,
    wrapper:
        "v6.0.0/bio/multiqc"


rule preprocess__multiqc__all:
    input:
        rules.preprocess__multiqc.output,
