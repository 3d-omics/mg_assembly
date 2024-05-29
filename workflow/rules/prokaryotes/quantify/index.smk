rule prokaryotes__quantify__index__:
    """Index dereplicader"""
    input:
        contigs=DREP / "dereplicated_genomes.fa.gz",
    output:
        mock=touch(QUANT_INDEX / "dereplicated_genomes"),
    log:
        QUANT_INDEX / "dereplicated_genomes.log",
    conda:
        "__environment__.yml"
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule prokaryotes__quantify__index:
    input:
        rules.prokaryotes__quantify__index__.output,
