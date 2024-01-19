rule _assemble__bowtie2__build:
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
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule assemble__bowtie2__build:
    """Index all megahit assemblies"""
    input:
        [ASSEMBLE_INDEX / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule _assemble__bowtie2__map:
    """Map one sample to one megahit assembly"""
    input:
        mock=ASSEMBLE_INDEX / "{assembly_id}",
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        reference=MEGAHIT / "{assembly_id}.fa.gz",
        fai=MEGAHIT / "{assembly_id}.fa.gz.fai",
    output:
        cram=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
    log:
        log=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
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
        | samtools sort \
            -l 9 \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2>> {log} 1>&2
        """


# rule _assemble__bowtie2__cram_to_bam:
#     input:
#         cram=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
#         reference=ASSEMBLE_RENAME / "{assembly_id}.fa",
#         fai=ASSEMBLE_RENAME / "{assembly_id}.fa.fai",
#     output:
#         bam=temp(ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam"),
#     log:
#         ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram_to_bam.log",
#     conda:
#         "__environment__.yml"
#     threads: 1
#     resources:
#         mem_mb=16 * 1024,
#         runtime=24 * 60,
#     shell:
#         """
#         samtools view \
#             --exclude-flags 4 \
#             --fast \
#             --output {output.bam} \
#             --output-fmt BAM \
#             --reference {input.reference} \
#             --threads {threads} \
#             {input.cram} \
#         2> {log} 1>&2
#         """


rule assemble__bowtie2:
    """Map all samples to all the assemblies that they belong to"""
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
