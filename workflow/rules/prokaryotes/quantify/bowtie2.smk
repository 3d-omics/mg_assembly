rule prokaryotes__quantify__bowtie2:
    """Align one sample to the dereplicated genomes"""
    input:
        mock=QUANT_INDEX / "drep.{secondary_ani}",
        forward_=PRE_BOWTIE2 / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "{sample_id}.{library_id}_2.fq.gz",
        reference=PROK_ANN / "drep.{secondary_ani}.fa.gz",
        fai=PROK_ANN / "drep.{secondary_ani}.fa.gz.fai",
    output:
        bam=QUANT_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.bam",
    log:
        QUANT_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        samtools_mem=params["quantify"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    shell:
        """
        find \
            $(dirname {output.bam}) \
            -name "$(basename {output.bam}).tmp.*.bam" \
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
            -M \
            -m {params.samtools_mem} \
            -o {output.bam} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2>> {log} 1>&2
        """


rule prokaryotes__quantify__bowtie2__all:
    """Align all samples to the dereplicated genomes"""
    input:
        [
            QUANT_BOWTIE2 / f"{secondary_ani}" / f"{sample_id}.{library_id}.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
            for secondary_ani in SECONDARY_ANIS
        ],
