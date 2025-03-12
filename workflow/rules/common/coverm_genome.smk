rule coverm__genome:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        "sample.{method}.bam",
    output:
        "sample.{method}.tsv.gz",  # it can be tsv too
    log:
        "sample.{method}.log",
    conda:
        "../../environments/coverm.yml"
    params:
        method=lambda w: w.method,
        separator="@",
        extra="--min-covered-fraction 0",
    resources:
        runtime: 60,
        mem_mb: 8 * 1024,
    shell:
        """
        ( coverm genome \
            --bam-files {input} \
            --methods {params.method} \
            --separator {params.separator} \
            {params.extra} \
        | cut \
            --fields 1 \
            --delimiter " " \
        | sed \
            '1 s/^Genome/sequence_id/g' \
        | gzip \
            --best \
        > {output} \
        ) 2> {log}
        """
