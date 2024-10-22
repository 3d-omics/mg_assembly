rule viruses__annotate__virsorter2:
    input:
        fna=GENOMADA / "rep_seq_virus.fna.gz",
        database=features["databases"]["virsorter2"],
    output:
        viruses_boundary=VIRSORTER2 / "final-viral-boundary.tsv.gz",
        combined=VIRSORTER2 / "final-viral-combined.fa.gz",
        score=VIRSORTER2 / "final-viral-score.tsv.gz",
        fa=VIRSORTER2 / "final-viral-combined-for-dramv.fa.gz",
        tsv=VIRSORTER2 / "viral-affi-contigs-for-dramv.tab.gz",
    log:
        VIRSORTER2 / "virsorter2.log",
    conda:
        "../../../environments/virsorter2.yml"
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
            {params.workdir}/for-dramv/viral-affi-contigs-for-dramv.tab \
            {params.workdir}/for-dramv/final-viral-combined-for-dramv.fa \
            {VIRSORTER2}/ \
        2>> {log} 1>&2

        bgzip \
            --threads {threads} \
            {VIRSORTER2}/final-viral-boundary.tsv \
            {VIRSORTER2}/final-viral-combined.fa \
            {VIRSORTER2}/final-viral-score.tsv \
            {VIRSORTER2}/final-viral-combined-for-dramv.fa \
            {VIRSORTER2}/viral-affi-contigs-for-dramv.tab \
        2>> {log} 1>&2
        """


rule viruses__annotate__virsorter2__all:
    input:
        rules.viruses__annotate__virsorter2.output,
