
rule viruses__quantify__bowtie2__:
    """Align one sample to the dereplicated genomes"""
    input:
        mock=VINDEX / "viruses",
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        reference=MMSEQS / "rep_seq.fasta.gz",
        fai=MMSEQS / "rep_seq.fasta.gz.fai",
    output:
        cram=VBOWTIE2 / "{sample_id}.{library_id}.cram",
    log:
        VBOWTIE2 / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        samtools_mem=params["quantify"]["samtools"]["mem"],
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


rule viruses__quantify__bowtie2:
    """Align all samples to the dereplicated genomes"""
    input:
        [
            VBOWTIE2 / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
