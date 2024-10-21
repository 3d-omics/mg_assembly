rule prokaryotes__quantify__index:
    """Index dereplicader"""
    input:
        contigs=PROK_ANN / "drep.{secondary_ani}.fa.gz",
    output:
        mock=touch(QUANT_INDEX / "drep.{secondary_ani}"),
    log:
        QUANT_INDEX / "drep.{secondary_ani}.log",
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


rule prokaryotes__quantify__index__all:
    input:
        [
            QUANT_INDEX / f"drep.{secondary_ani}.fa.gz"
            for secondary_ani in SECONDARY_ANIS
        ],
