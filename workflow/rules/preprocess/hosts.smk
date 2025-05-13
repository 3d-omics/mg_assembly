rule preprocess__hosts:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=lambda wildcards: features["hosts"][wildcards.host],
    output:
        fa_gz=PRE_HOSTS / "{host}.fa.gz",
    log:
        PRE_HOSTS / "{host}.log",
    conda:
        "../../environments/hosts.yml"
    cache: "omit-software"
    threads: 8
    shell:
        """
        ( gzip \
            --decompress \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule preprocess__hosts__all:
    """Recompress all host genomes"""
    input:
        [PRE_HOSTS / f"{host}.fa.gz" for host in HOST_NAMES],
