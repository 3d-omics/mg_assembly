rule assemble_bowtie2_build_one:
    """Index a megahit assembly"""
    input:
        contigs=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        mock=touch(ASSEMBLE_INDEX / "{assembly_id}"),
    log:
        ASSEMBLE_INDEX / "{assembly_id}.log",
    conda:
        "assemble.yml"
    threads: 24
    params:
        extra=params["assemble"]["bowtie2-build"]["extra"],
    resources:
        mem_mb=double_ram(params["assemble"]["bowtie2-build"]["memory_gb"]),
        runtime=48 * 60,
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule assemble_bowtie2_build_all:
    """Index all megahit assemblies"""
    input:
        [ASSEMBLE_INDEX / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule assemble_bowtie2_one:
    """Map one sample to one megahit assembly"""
    input:
        mock=ASSEMBLE_INDEX / "{assembly_id}",
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        reference=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        cram=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
    log:
        log=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        "assemble.yml"
    threads: 24
    params:
        extra=params["assemble"]["bowtie2"]["extra"],
        samtools_mem=params["assemble"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    resources:
        mem_mb=double_ram(params["assemble"]["bowtie2"]["memory_gb"]),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        find \
            $(dirname {output.cram}) \
            -name "$(basename {output.cram}).tmp.*.bam" \
            -delete \
        2> {log} 1>&2

        ( bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2>> {log} 1>&2
        """


rule assemble_bowtie2:
    """Map all samples to all the assemblies that they belong to"""
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
