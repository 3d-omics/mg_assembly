rule coverm__contig:
    """Compute coverm statistics for a single bam using the method {method}."""
    input:
        "sample.{method}.bam",
    output:
        "sample.{method}.tsv.gz",
    log:
        "sample.{method}.log",
    conda:
        "../../environments/coverm.yml"
    params:
        method=lambda w: w.method,
    shell:
        """
        ( coverm contig \
            --bam-files {input} \
            --methods {params.method} \
            --proper-pairs-only \
        | cut \
            --fields 1 \
            --delimiter " " \
        | sed \
            '1 s/^Contig/sequence_id/g' \
        | gzip \
            --best \
        > {output} \
        ) 2> {log}
        """
