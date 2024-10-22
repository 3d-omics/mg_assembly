rule viruses__multiqc:
    input:
        bowtie2=[
            VBOWTIE2 / f"{sample_id}.{library_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in BAM_REPORTS
        ],
        quast=rules.viruses__annotate__quast__all.input,
    output:
        html=RESULTS / "viruses.html",
        folder=directory(RESULTS / "viruses_data"),
    log:
        RESULTS / "viruses.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        outdir=RESULTS,
    shell:
        """
        multiqc \
            --title viruses \
            --force \
            --filename viruses \
            --outdir {params.outdir} \
            --dirs \
            --dirs-depth 1 \
            --fullnames \
            {input} \
        2> {log} 1>&2
        """


rule viruses__multiqc__all:
    input:
        rules.viruses__multiqc.output,
