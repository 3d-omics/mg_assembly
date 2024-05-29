rule prokaryotes__annotate__quast__:
    """Run quast over one the dereplicated mags"""
    input:
        DREP / "dereplicated_genomes.fa.gz",
    output:
        directory(QUAST),
    log:
        QUAST / "quast.log",
    conda:
        "__environment__.yml"
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """


rule prokaryotes__annotate__quast:
    input:
        rules.prokaryotes__annotate__quast__.output,
