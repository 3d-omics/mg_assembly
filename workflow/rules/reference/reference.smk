rule reference_recompress:
    input:
        fa_gz=lambda wildcards: features["hosts"][wildcards.genome],
    output:
        REFERENCE / "{genome}.fa.gz",
    log:
        REFERENCE / "{genome}.log",
    threads: 24
    conda:
        "reference.yml"
    shell:
        """
        (gzip -dc {input.fa_gz} | bgzip -@ {threads} > {output}) 2> {log}
        """


rule reference:
    input:
        [REFERENCE / f"{genome}.fa.gz" for genome in HOST_NAMES],
