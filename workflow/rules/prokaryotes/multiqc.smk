rule prokaryotes__multiqc:
    input:
        bowtie2=[
            PROK_BOWTIE2
            / f"drep.{secondary_ani}"
            / f"{sample_id}.{library_id}.stats.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
            for secondary_ani in SECONDARY_ANIS
        ],
        quast=[
            PROK_QUAST / f"drep.{secondary_ani}"
            for secondary_ani in SECONDARY_ANIS
        ],
    output:
        RESULTS / "prokaryotes.html",
        RESULTS / "prokaryotes_data.zip",
    log:
        RESULTS / "prokaryotes.log",
    params:
        extra="--title prokaryotes --dirs --fullnames --fn_as_s_name --force",
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=6 * 60,
    wrapper:
        "v6.0.0/bio/multiqc"


rule prokaryotes__multiqc__all:
    input:
        rules.prokaryotes__multiqc.output,
