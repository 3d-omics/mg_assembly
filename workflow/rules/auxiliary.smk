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
