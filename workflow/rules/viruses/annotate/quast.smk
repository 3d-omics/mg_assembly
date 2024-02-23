rule viruses__annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        MMSEQS / "rep_seq.fasta",
    output:
        directory(QUASTV),
    log:
        QUASTV / "quast.log",
    conda:
        "__environment__.yml"
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=1 * 60,
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """
