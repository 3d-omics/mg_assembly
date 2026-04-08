rule concatenate_flat_ungzipped:
    """Concatenate multiple flat files (fasta, txt, gff, genbank) and recompress into a single one"""
    input:
        ["fasta1.fa", "fasta2.fa"],
    output:
        "fasta_out.fa.gz",
    conda:
        ENVS / "concatenate.yml"
    threads: 24
    log:
        "fasta_out.log",
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


rule concatenate_tsv_ungzipped:
    """
    Concatenate multiple tsvs respecting column order, and compress into a single file.
    """
    input:
        ["table1.tsv", "table2.tsv"],
    output:
        "table.tsv.gz",
    threads: 24
    log:
        "table.log",
    params:
        compress_level=5,
    conda:
        ENVS / "concatenate.yml"
    shell:
        """
        (
            csvtk concat \
                --tabs \
                {input} /dev/null \
            | bgzip \
                --compress-level {params.compress_level} \
                --threads {threads} \
            > {output}
        ) 2> {log}
        """
