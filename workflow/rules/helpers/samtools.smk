rule samtools_index_crai:
    """Index a cram file"""
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    conda:
        "_env.yml"
    log:
        "{prefix}.cram.crai.log",
    shell:
        "samtools index {input} 2> {log}"


rule samtools_index_bai:
    """Index a bam file"""
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    conda:
        "_env.yml"
    log:
        "{prefix}.bam.bai.log",
    shell:
        "samtools index {input} 2> {log}"


rule samtools_faidx_fa:
    """Index a fasta file"""
    input:
        "{prefix}.fa",
    output:
        "{prefix}.fa.fai",
    conda:
        "_env.yml"
    log:
        "{prefix}.fa.fai.log",
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule samtools_faidx_fa_gz:
    """Index a gzipped fasta file"""
    input:
        "{prefix}.fa.gz",
    output:
        fai="{prefix}.fa.gz.fai",
        gzi="{prefix}.fa.gz.gzi",
    conda:
        "_env.yml"
    log:
        "{prefix}.fa.gz.log",
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule samtools_idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "_env.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"


rule samtools_flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "_env.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"
