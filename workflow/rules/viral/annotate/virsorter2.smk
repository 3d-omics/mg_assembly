rule _viral__annotate__virsorter2:
    input:
        fna=GENOMADA / "rep_seq_virus.fna",
        database=features["databases"]["virsorter2"],
    output:
        viral_boundary=VIRSORTER2 / "final-viral-boundary.tsv",
        combined=VIRSORTER2 / "final-viral-combined.fa",
        score=VIRSORTER2 / "final-viral-score.tsv",
        fa=VIRSORTER2 / "final-viral-combined-for-dramv.fa",
        tsv=VIRSORTER2 / "viral-affi-contigs-for-dramv.tab",
    log:
        VIRSORTER2 / "virsorter2.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        workdir=VIRSORTER2,
    # shadow: "minimal"
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


rule viral__annotate__virsorter2:
    input:
        VIRSORTER2 / "for-dramv" / "final-viral-combined-for-dramv.fa",
