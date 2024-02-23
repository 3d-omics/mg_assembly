rule prokaryotes__cluster__index__:
    """Index a megahit assembly"""
    input:
        contigs=MEGAHIT / "{assembly_id}.fa.gz",
    output:
        mock=touch(ASSEMBLE_INDEX / "{assembly_id}"),
    log:
        ASSEMBLE_INDEX / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["assemble"]["bowtie2-build"]["memory_gb"]),
        runtime=48 * 60,
        attempt=get_attempt,
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.contigs} \
            {output.mock} \
        2> {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """


rule prokaryotes__cluster__index:
    """Index all megahit assemblies"""
    input:
        [ASSEMBLE_INDEX / f"{assembly_id}" for assembly_id in ASSEMBLIES],

