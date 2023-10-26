rule assemble_samtools_stats_cram_assembly:
    """Run samtools stats over one library, using as reference an assembly"""
    input:
        cram=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
        crai=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram.crai",
        reference=ASSEMBLE_RENAME / "{assembly_id}.fa",
        fai=ASSEMBLE_RENAME / "{assembly_id}.fa.fai",
    output:
        txt=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.stats.txt",
    log:
        ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.stats.log",
    conda:
        "assemble.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


rule assemble_samtools:
    """Get all the samtools stats for all assemblies and samples"""
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.{extension}"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
        ],
