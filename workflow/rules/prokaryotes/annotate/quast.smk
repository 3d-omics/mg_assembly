rule prokaryotes__annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        PROK_ANN / "drep.{secondary_ani}.fa.gz",
    output:
        directory(QUAST / "drep.{secondary_ani}"),
    log:
        QUAST / "drep.{secondary_ani}.log",
    conda:
        "../../../environments/quast.yml"
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
        [
            PROK_ANN / f"drep.{secondary_ani}.fa.gz"
            for secondary_ani in SECONDARY_ANIS
        ],
