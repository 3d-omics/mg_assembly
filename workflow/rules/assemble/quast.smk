rule assemble_eval_quast_one:
    """Run quast over one assembly group"""
    input:
        ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        directory(ASSEMBLE_QUAST / "{assembly_id}"),
    log:
        ASSEMBLE_QUAST / "{assembly_id}.log",
    conda:
        "assemble.yml"
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


rule assemble_eval_quast:
    """Run quast over all assembly groups"""
    input:
        [ASSEMBLE_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],