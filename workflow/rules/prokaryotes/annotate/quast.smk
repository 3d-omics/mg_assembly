rule prokaryotes__annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        [PROK_ANN / f"drep.{secondary_ani}.fa.gz" for secondary_ani in SECONDARY_ANIS],
    output:
        directory(PROK_QUAST),
    log:
        PROK_ANN / "quast.log",
    conda:
        "../../../environments/quast.yml"
    threads: 4
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """


rule prokaryotes__annotate__quast__all:
    input:
        [PROK_ANN / f"drep.{secondary_ani}.fa.gz" for secondary_ani in SECONDARY_ANIS],
