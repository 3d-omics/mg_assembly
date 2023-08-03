include: "helpers_functions.smk"


rule crai:
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    conda:
        "../envs/auxiliary.yml"
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
        "../envs/auxiliary.yml"
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
        "../envs/auxiliary.yml"
    log:
        "{prefix}.fa.fai.log",
    shell:
        "samtools faidx {input} 2> {log} 1>&2"
