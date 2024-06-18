rule prokaryotes__quantify__bowtie2__:
    """Align one sample to the dereplicated genomes"""
    input:
        mock=QUANT_INDEX / "dereplicated_genomes",
        forward_=get_host_clean_fastq_forward,
        reverse_=get_host_clean_fastq_reverse,
        reference=DREP / "dereplicated_genomes.fa.gz",
        fai=DREP / "dereplicated_genomes.fa.gz.fai",
    output:
        cram=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram",
    log:
        QUANT_BOWTIE2 / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        samtools_mem=params["quantify"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
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
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2>> {log} 1>&2
        """


rule prokaryotes__quantify__bowtie2:
    """Align all samples to the dereplicated genomes"""
    input:
        [
            QUANT_BOWTIE2 / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
