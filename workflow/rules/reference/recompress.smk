rule reference__recompress__:
    """Recompress the reference with bgzip"""
    input:
        fa_gz=lambda wildcards: features["hosts"][wildcards.genome],
    output:
        HOSTS / "{genome}.fa.gz",
    log:
        HOSTS / "{genome}.log",
    conda:
        "__environment__.yml"
    cache: True
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.fa_gz} \
        | bgzip \
            --threads {threads} \
        > {output} \
        ) 2> {log}
        """


rule reference__recompress:
    """Run all the steps for reference(s) preparations"""
    input:
        [HOSTS / f"{genome}.fa.gz" for genome in HOST_NAMES],
