rule fastqc:
    input:
        "{prefix}.fq.gz",
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip",
    conda:
        "helpers.yml"
    log:
        "{prefix}_fastqc.log",
    shell:
        "fastqc {input} 2> {log} 1>&2"


rule crai:
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    conda:
        "helpers.yml"
    log:
        "{prefix}.cram.crai.log",
    shell:
        "samtools index {input} 2> {log}"


rule bai:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    conda:
        "helpers.yml"
    log:
        "{prefix}.bam.bai.log",
    shell:
        "samtools index {input} 2> {log}"


rule fa_fai:
    input:
        "{prefix}.fa",
    output:
        "{prefix}.fa.fai",
    conda:
        "helpers.yml"
    log:
        "{prefix}.fa.fai.log",
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
        "helpers.yml"
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
        "helpers.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"
