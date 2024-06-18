rule preprocess__samtools__stats_cram_host:
    """Compute stats for a host cram"""
    input:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
        crai=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram.crai",
        reference=HOSTS / "{genome}.fa.gz",
    output:
        tsv=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.stats.tsv",
    log:
        PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output.tsv} \
        2> {log}
        """
