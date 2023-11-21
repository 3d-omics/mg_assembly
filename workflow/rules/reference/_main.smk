rule _reference__recompress:
    """Recompress the reference with bgzip"""
    input:
        fa_gz=lambda wildcards: features["hosts"][wildcards.genome],
    output:
        REFERENCE / "{genome}.fa.gz",
    log:
        REFERENCE / "{genome}.log",
    threads: 24
    conda:
        "_env.yml"
    shell:
        """
        (gzip -dc {input.fa_gz} | bgzip -@ {threads} > {output}) 2> {log}
        """


rule reference:
    """Run all the steps for reference(s) preparations"""
    input:
        [REFERENCE / f"{genome}.fa.gz" for genome in HOST_NAMES],
