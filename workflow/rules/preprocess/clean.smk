rule preprocess__clean__cram_to_fastq:
    input:
        cram=get_host_clean_cram,
    output:
        forward_=CLEAN / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=CLEAN / "{sample_id}.{library_id}_2.fq.gz",
    log:
        CLEAN / "{sample_id}.{library_id}.cram_to_fastq.log",
    conda:
        "__environment__.yml"
    shell:
        """
        rm -rf {output.forward_}.collate

        ( samtools view \
            -f 12 \
            -u \
            --threads {threads} \
            {input} \
        | samtools collate \
            -O \
            -u \
            -T {output.forward_}.collate \
            --threads {threads} \
            - \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            --threads {threads} \
            -c 0 \
            /dev/stdin \
        ) 2> {log} 1>&2
        """


rule preprocess__clean:
    input:
        [
            CLEAN / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
