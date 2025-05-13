rule assemble__multiqc:
    input:
        bowtie2=[
            ASMB_BOWTIE2 / assembly_id / f"{sample_id}.{library_id}.{report}"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
            for report in BAM_REPORTS
        ],
        quast=[ASMB_QUAST / assembly_id for assembly_id in ASSEMBLIES],
        kraken2=[
            ASMB_KRAKEN2 / kraken2_db / f"{assembly_id}.report"
            for kraken2_db in features["databases"]["kraken2"]
            for assembly_id in ASSEMBLIES
        ],
    output:
        RESULTS / "assemble.html",
        RESULTS / "assemble_data.zip",
    log:
        RESULTS / "assemble.log",
    params:
        extra="--title assemble --dirs --fullnames --fn_as_s_name --force",
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=6 * 60,
    wrapper:
        "v6.0.0/bio/multiqc"


rule assemble__multiqc__all:
    input:
        rules.assemble__multiqc.output,
