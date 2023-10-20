rule dereplicate_samtools_stats_cram_dereplicate:
    input:
        cram=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa",
    output:
        txt=DREP_BOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        DREP_BOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "helpers.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


rule dereplicate_eval_samtools:
    input:
        [
            DREP_BOWTIE2 / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
        ],
