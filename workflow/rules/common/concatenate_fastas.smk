rule concatenate_fastas:
    input:
        ["fasta1.fa.gz", "fasta2.fa.gz"]
    output:
        "fasta_out.fa.gz"
    conda:
        "../../environments/htslib.yml"
    threads:
        1
    log:
        "fasta_out.log"
    params:
        compress_level = 5,
    shell:
        """
        ( bgzip \
            --decompress \
            --stdout \
            {input} \
        | bgzip \
            --compress-level {params.compress_level} \
            --threads {threads} \
        > {output} \
        ) 2> {log}
        """
    