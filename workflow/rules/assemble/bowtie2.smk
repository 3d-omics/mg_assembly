rule assemble__bowtie2:
    """Map one sample to one megahit assembly"""
    input:
        mock=ASSEMBLE_INDEX / "{assembly_id}",
        forward_=PRE_BOWTIE2 / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "{sample_id}.{library_id}_2.fq.gz",
        reference=ASSEMBLE_MEGAHIT / "{assembly_id}.fa.gz",
        fai=ASSEMBLE_MEGAHIT / "{assembly_id}.fa.gz.fai",
    output:
        bam=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
    log:
        log=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        samtools_mem=params["assemble"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    resources:
        attempt=get_attempt,
    retries: 5
    shell:
        """
        find \
            $(dirname {output.bam}) \
            -name "$(basename {output.bam}).tmp.*.bam" \
            -delete \
        2> {log}.{resources.attempt} 1>&2

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
            -o {output.bam} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2>> {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """


rule assemble__bowtie2__all:
    """Map all samples to all the assemblies that they belong to"""
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ]
        + [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam.bai"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
