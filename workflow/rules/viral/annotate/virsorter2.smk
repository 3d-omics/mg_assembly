rule _viral__annotate__virsorter2:
    input:
        fna=GENOMADA / "genomad_virus.fna",
        database=features["databases"]["virsorter2"],
    output:
        fa=VIRSORTER2 / "for-dramv" / "final-viral-combined-for-dramv.fa",
        tsv=VIRSORTER2 / "for-dramv" / "viral-affi-contigs-for-dramv.tab",
    log:
        VIRSORTER2 / "virsorter2.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        workdir=VIRSORTER2,
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
        """


rule viral__annotate__virsorter2:
    input:
        VIRSORTER2 / "for-dramv" / "final-viral-combined-for-dramv.fa",
