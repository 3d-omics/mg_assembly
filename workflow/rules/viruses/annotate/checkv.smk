rule viruses__annotate__checkv:
    input:
        fasta=MMSEQS / "rep_seq.fa.gz",
        database=features["databases"]["checkv"],
    output:
        complete_genomes=CHECKV / "complete_genomes.tsv",
        completeness=CHECKV / "completeness.tsv",
        contamination=CHECKV / "contamination.tsv",
        proviruses=CHECKV / "proviruses.fna",
        summary=CHECKV / "quality_summary.tsv",
        viruses=CHECKV / "viruses.fna",
    log:
        CHECKV / "checkv.log",
    conda:
        "../../../environments/checkv.yml"
    params:
        workdir=CHECKV,
    shell:
        """
        checkv end_to_end \
            -d {input.database} \
            -t {threads} \
            --restart \
            {input.fasta} \
            {params.workdir} \
        2> {log} 1>&2
        """


rule viruses__annotate__checkv__all:
    input:
        rules.viruses__annotate__checkv.output,
