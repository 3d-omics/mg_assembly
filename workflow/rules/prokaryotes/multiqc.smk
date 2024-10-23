rule prokaryotes__multiqc:
    input:
        bowtie2=[
            QUANT_BOWTIE2
            / f"drep.{secondary_ani}"
            / f"{sample_id}.{library_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in ["stats.tsv", "flagstats.txt"]
            for secondary_ani in SECONDARY_ANIS
        ],
        quast=rules.prokaryotes__annotate__quast__all.output,
    output:
        html=RESULTS / "prokaryotes.html",
        folder=directory(RESULTS / "prokaryotes_data"),
    log:
        RESULTS / "prokaryotes.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        outdir=RESULTS,
    shell:
        """
        multiqc \
            --title assemble \
            --force \
            --filename assemble \
            --outdir {params.outdir} \
            --dirs \
            --dirs-depth 1 \
            --fullnames \
            {input} \
        2> {log} 1>&2
        """


rule prokaryotes__multiqc__all:
    input:
        rules.prokaryotes__multiqc.output,
