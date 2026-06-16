rule viruses__multiqc:
    input:
        bowtie2=[
            VIR_BOWTIE2 / "rep_seq" / f"{sample_id}.{library_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in BAM_REPORTS
        ],
        quast=VIR_QUAST,
    output:
        RESULTS / "viruses.html",
        RESULTS / "viruses_data.zip",
    log:
        RESULTS / "viruses.log",
    params:
        extra="--title viruses --dirs --fullnames --fn_as_s_name --force",
    wrapper:
        "v9.4.0/bio/multiqc"


rule viruses__multiqc__all:
    input:
        rules.viruses__multiqc.output,
