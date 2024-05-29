rule viruses__annotate__virsorter2__:
    input:
        fna=GENOMADA / "rep_seq_virus.fna",
        database=features["databases"]["virsorter2"],
    output:
        viruses_boundary=VIRSORTER2 / "final-viruses-boundary.tsv",
        combined=VIRSORTER2 / "final-viruses-combined.fa",
        score=VIRSORTER2 / "final-viruses-score.tsv",
        fa=VIRSORTER2 / "final-viruses-combined-for-dramv.fa",
        tsv=VIRSORTER2 / "viruses-affi-contigs-for-dramv.tab",
    log:
        VIRSORTER2 / "virsorter2.log",
    conda:
        "__environment__.yml"
    params:
        workdir=VIRSORTER2,
    shadow:
        "minimal"
    shell:
        """
        virsorter run \
            --working-dir {params.workdir} \
            --jobs {threads} \
            --prep-for-dramv \
            --tmpdir {params.workdir}/tmp \
            --rm-tmpdir \
            --verbose \
            --seqfile {input.fna} \
            --db-dir {input.database} \
        2> {log} 1>&2

        mv \
            {params.workdir}/for-dramv/* \
            {params.workdir}/
        """


rule viruses__annotate__virsorter2:
    input:
        rules.viruses__annotate__virsorter2__.output,
