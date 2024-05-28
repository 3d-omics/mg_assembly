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
    shell:
        """
        (gzip -dc {input.fa_gz} | bgzip -@ {threads} > {output}) 2> {log}
        """


rule reference__recompress:
    """Run all the steps for reference(s) preparations"""
    input:
        [HOSTS / f"{genome}.fa.gz" for genome in HOST_NAMES],
