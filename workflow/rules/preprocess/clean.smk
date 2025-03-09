include: "clean_functions.smk"


rule preprocess__clean:
    """Get the final fastq files and compress it properly"""
    input:
        forward_=get_final_fastq_forward,
        reverse_=get_final_fastq_reverse,
    output:
        forward_=PRE_CLEAN / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_CLEAN / "{sample_id}.{library_id}_2.fq.gz",
    log:
        PRE_CLEAN / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/bowtie2.yml"  # It has htslib in it
    threads: 24
    resources:
        mem_mb=1 * 1024,
        runtime=1 * 60,
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.forward_} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.forward_} \
        ) 2> {log}

        ( gzip \
            --decompress \
            --stdout \
            {input.reverse_} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.reverse_} \
        ) 2>> {log}
        """


rule preprocess__clean__all:
    input:
        [
            PRE_CLEAN / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
