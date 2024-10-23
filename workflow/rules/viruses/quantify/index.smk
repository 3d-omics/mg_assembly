rule viruses__quantify__index:
    """Index dereplicader"""
    input:
        contigs=MMSEQS / "rep_seq.fa.gz",
    output:
        mock=touch(VINDEX / "viruses"),
    log:
        VINDEX / "virues.log",
    conda:
        "../../../environments/bowtie2_samtools.yml"
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule viruses__quantify__index__all:
    input:
        rules.viruses__quantify__index.output,
