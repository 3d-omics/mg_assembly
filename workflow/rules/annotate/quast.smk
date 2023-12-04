rule annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        DREP / "dereplicated_genomes.fa",
    output:
        directory(QUAST),
    log:
        QUAST / "quast.log",
    conda:
        "_env.yml"
    threads: 4
    params:
        extra=params["dereplicate"]["quast"]["extra"],
    resources:
        mem_mb=8 * 1024,
        runtime=1 * 60,
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """
