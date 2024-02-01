rule _viral__annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        MMSEQS / "cluster.fa",
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


rule viral__annotate__quast:
    input:
        rules._viral__annotate__quast.output,
