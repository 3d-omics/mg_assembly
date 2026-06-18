rule quast:
    """Run quast over an assembly"""
    input:
        "assembly.fa.gz",
    output:
        directory("quast"),
    log:
        "quast.log",
    conda:
        ENVS / "quast.yml"
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=6 * 60,
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """
