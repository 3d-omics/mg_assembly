rule pre_samtools_stats_cram:
    """Compute the stats of a cram file using samtools stats"""
    input:
        cram=PRE_BOWTIE2 / "{genome}/{sample_id}.{library_id}.cram",
        crai=PRE_BOWTIE2 / "{genome}/{sample_id}.{library_id}.cram.crai",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        txt=PRE_BOWTIE2 / "{genome}/{sample_id}.{library_id}.stats.txt",
    log:
        PRE_BOWTIE2 / "{genome}/{sample_id}.{library_id}.stats.log",
    conda:
        "pre.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


rule pre_samtools:
    """Get all the stats of a bam file using samtools"""
    input:
        [
            PRE_BOWTIE2 / f"{genome}/{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
            for genome in HOST_NAMES
        ],
