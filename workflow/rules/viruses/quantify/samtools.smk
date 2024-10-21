rule viruses__quantify__samtools__stats_bam__:
    """Get stats from bam files using samtools stats."""
    input:
        bam=VBOWTIE2 / "{sample_id}.{library_id}.bam",
        bai=VBOWTIE2 / "{sample_id}.{library_id}.bam.bai",
        reference=MMSEQS / "rep_seq.fa.gz",
        fai=MMSEQS / "rep_seq.fa.gz.fai",
    output:
        txt=VBOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        VBOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools stats --reference {input.reference} {input.bam} > {output.txt} 2> {log}"


rule viruses__quantify__samtools:
    """Get stats from bam files using samtools stats."""
    input:
        [
            VBOWTIE2 / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt"]
        ],
