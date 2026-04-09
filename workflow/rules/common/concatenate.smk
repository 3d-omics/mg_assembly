rule concatenate__flat_to_gzipped:
    """Concatenate multiple flat files (fasta, txt, gff, genbank) and recompress into a single one"""
    input:
        ["fasta1.fa", "fasta2.fa"],
    output:
        "fasta_out.fa.gz",
    log:
        "fasta_out.log",
    conda:
        ENVS / "concatenate.yml"
    threads: 24
    params:
        compress_level=5,
    shell:
        """
        ( cat \
            {input} \
        | bgzip \
            --compress-level {params.compress_level} \
            --threads {threads} \
        > {output} \
        ) 2> {log}
        """


rule concatenate__gzipped_to_gzipped:
    """Concatenate multiple gzipped files (fasta, txt, gff, genbank) and recompress into a single one"""
    input:
        ["fasta1.fa.gz", "fasta2.fa.gz"],
    output:
        "fasta_out.fa.gz",
    log:
        "fasta_out.log",
    conda:
        ENVS / "concatenate.yml"
    threads: 24
    params:
        compress_level=5,
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input} \
        | bgzip \
            --compress-level {params.compress_level} \
            --threads {threads} \
        > {output} \
        ) 2> {log}
        """
