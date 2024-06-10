rule viruses__quantify__samtools__stats_cram__:
    """Get stats from CRAM files using samtools stats."""
    input:
        cram=VBOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=VBOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=MMSEQS / "rep_seq.fasta.gz",
        fai=MMSEQS / "rep_seq.fasta.gz.fai",
    output:
        txt=VBOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        VBOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


rule viruses__quantify__samtools:
    """Get stats from CRAM files using samtools stats."""
    input:
        [
            VBOWTIE2 / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt"]
        ],
