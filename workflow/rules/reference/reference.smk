rule reference_recompress:
    input:
        fa_gz=features["host"]["fasta"],
    output:
        REFERENCE / "host.fa.gz",
    log:
        REFERENCE / "host.log",
    threads: 24
    conda:
        "reference.yml"
    shell:
        """
        (gzip -dc {input.fa_gz} | bgzip -@ {threads} > {output}) 2> {log}
        """
