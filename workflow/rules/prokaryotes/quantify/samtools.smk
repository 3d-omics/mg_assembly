rule prokaryotes__quantify__samtools__stats_cram__:
    """Get stats from CRAM files using samtools stats."""
    input:
        cram=QUANT_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.cram",
        crai=QUANT_BOWTIE2 / "drep.{secondary_ani}" /"{sample_id}.{library_id}.cram.crai",
        reference=PROK_ANN / "drep.{secondary_ani}.fa.gz",
        fai=PROK_ANN / "drep.{secondary_ani}.fa.gz.fai",
    output:
        txt=QUANT_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.stats.txt",
    log:
        QUANT_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output.txt} \
        2> {log}
        """


rule prokaryotes__quantify__samtools:
    """Get stats from CRAM files using samtools stats."""
    input:
        [
            QUANT_BOWTIE2 / f"drep.{secondary_ani}" / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt"]
            for secondary_ani in SECONDARY_ANIS
        ],
