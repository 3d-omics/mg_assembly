use rule multiqc as viruses__multiqc with:
    input:
        bowtie2=[
            VIR_MVP_READ_MAPPING / f"{sample_id}" / f"{sample_id}_sorted.stats.tsv"
            for sample_id in SAMPLES
        ],
        quast=VIR_QUAST,
    output:
        RESULTS / "viruses.html",
        RESULTS / "viruses_data.zip",
    log:
        RESULTS / "viruses.log",
    params:
        extra="--title viruses --dirs --fullnames --fn_as_s_name --force",


rule viruses__multiqc__all:
    input:
        rules.viruses__multiqc.output,
