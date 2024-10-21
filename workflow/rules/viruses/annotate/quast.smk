rule viruses__annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        MMSEQS / "rep_seq.fa.gz",
    output:
        directory(QUASTV),
    log:
        QUASTV / "quast.log",
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


rule viruses__annotate__quast__all:
    input:
        rules.viruses__annotate__quast.output,
