rule prokaryotes__annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        PROK_ANN / "drep.{secondary_ani}.fa.gz",
    output:
        directory(PROK_QUAST / "drep.{secondary_ani}"),
    log:
        PROK_QUAST / "drep.{secondary_ani}.log",
    conda:
        ENVS / "quast.yml"
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
        [PROK_QUAST / f"drep.{secondary_ani}" for secondary_ani in SECONDARY_ANIS],
