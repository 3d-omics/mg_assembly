rule dereplicate_quast:
    """Run quast over one the dereplicated mags"""
    input:
        DREP / "dereplicated_genomes.fa",
    output:
        directory(DREP_QUAST / "dereplicated_genomes"),
    log:
        DREP_QUAST / "dereplicated_genomes.log",
    conda:
        "dereplicate.yml"
    threads: 4
    params:
        extra=params["dereplicate"]["quast"]["extra"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """
