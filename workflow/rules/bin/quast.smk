rule _bin__quast__run:
    """Run quast over one assembly group"""
    input:
        MAGSCOT / "{assembly_id}.fa",
    output:
        directory(BIN_QUAST / "{assembly_id}"),
    log:
        BIN_QUAST / "{assembly_id}.log",
    conda:
        "_env.yml"
    threads: 4
    params:
        extra=params["assemble"]["quast"]["extra"],
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


rule bin__quast:
    """Run quast over all assembly groups"""
    input:
        [BIN_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],
