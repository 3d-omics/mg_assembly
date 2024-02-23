rule preprocess__index__:
    """Build PRE_BOWTIE2 index for the host reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=HOSTS / "{genome}.fa.gz",
        faidx=HOSTS / "{genome}.fa.gz.fai",
    output:
        mock=touch(PRE_INDEX / "{genome}"),
    log:
        PRE_INDEX / "{genome}.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["preprocess"]["bowtie2-build"]["memory_gb"]),
        runtime=24 * 60,
        attempt=get_attempt,
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.reference} \
            {output.mock} \
        2> {log}.{resources.attempt} 1>&2

        mv \
            {log}.{resources.attempt} \
            {log}
        """


rule preprocess__index:
    input:
        [PRE_INDEX / f"{genome}" for genome in HOST_NAMES],
