rule assemble__bowtie2__build:
    """Index a megahit assembly"""
    input:
        contigs=ASSEMBLE_MEGAHIT / "{assembly_id}.fa.gz",
    output:
        mock=touch(ASSEMBLE_INDEX / "{assembly_id}"),
    log:
        ASSEMBLE_INDEX / "{assembly_id}.log",
    conda:
        "../../environments/bowtie2_samtools.yml"
    resources:
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


rule assemble__bowtie2__build__all:
    """Index all megahit assemblies"""
    input:
        [ASSEMBLE_INDEX / f"{assembly_id}" for assembly_id in ASSEMBLIES],
