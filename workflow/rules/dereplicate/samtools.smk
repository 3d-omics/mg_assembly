rule _dereplicate__samtools__stats_cram:
    """Get stats from CRAM files using samtools stats."""
    input:
        cram=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa",
        fai=DREP / "dereplicated_genomes.fa.fai",
    output:
        txt=DREP_BOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        DREP_BOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "_env.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


rule dereplicate__samtools:
    """Get stats from CRAM files using samtools stats."""
    input:
        [
            DREP_BOWTIE2 / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
        ],
