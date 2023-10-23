rule bin_quast_one:
    """Run quast over one assembly group"""
    input:
        MAGSCOT / "{assembly_id}.fa",
    output:
        directory(BIN_QUAST / "{assembly_id}"),
    log:
        BIN_QUAST / "{assembly_id}.log",
    conda:
        "bin.yml"
    threads: 4
    params:
        extra=params["assemble"]["quast"]["extra"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """


rule bin_quast:
    """Run quast over all assembly groups"""
    input:
        [BIN_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],
