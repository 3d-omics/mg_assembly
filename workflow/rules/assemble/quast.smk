rule assemble__quast:
    """Run quast over each assembly"""
    input:
        ASMB_MEGAHIT / "{assembly_id}.fa.gz",
    output:
        directory(ASMB_QUAST / "{assembly_id}"),
    log:
        ASMB_QUAST / "{assembly_id}.log",
    conda:
        "../../environments/quast.yml"
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """


rule assemble__quast__all:
    input:
        [ASMB_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],
