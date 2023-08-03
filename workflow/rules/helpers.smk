include: "helpers_functions.smk"


rule fastqc:
    input:
        "{prefix}.fq.gz",
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip",
    conda:
        "../envs/helpers.yml"
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
        "../envs/helpers.yml"
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
        "../envs/helpers.yml"
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
        "../envs/helpers.yml"
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
        "../envs/helpers.yml"
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
        "../envs/helpers.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule samtools_stats_cram_pre:
    input:
        cram=PRE_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=PRE_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=features["host"]["fasta"],
    output:
        txt=PRE_BOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        PRE_BOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "../envs/helpers.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


rule samtools_stats_cram_assembly:
    input:
        cram=ASSEMBLY_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.cram",
        crai=ASSEMBLY_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.cram.crai",
        reference=ASSEMBLY_RENAME / "{assemble_id}.fa",
    output:
        txt=ASSEMBLY_BOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        ASSEMBLY_BOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "../envs/helpers.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


# rule samtools_stat_cram_bin


rule samtools_stats_cram_metabin:
    input:
        cram=METABIN_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.cram",
        crai=METABIN_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.cram.crai",
        reference=METABIN_RENAME / "{assemble_id}.fa",
    output:
        txt=METABIN_BOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        METABIN_BOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "../envs/helpers.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


rule samtools_stats_cram_dereplicate:
    input:
        cram=DREP_BOWTIE2 / "dereplicated_genomes.{sample_id}.{library_id}.cram",
        crai=DREP_BOWTIE2 / "dereplicated_genomes.{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa",
    output:
        txt=DREP_BOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        DREP_BOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "../envs/helpers.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"
