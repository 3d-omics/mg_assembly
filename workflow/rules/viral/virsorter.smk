rule _viral__virsorer__run:
    input:
        fna=CHECKV / "{assembly_id}" / "all.fna",
        database=features["databases"]["virsorter2"],
    output:
        fa=VIRSORTER2
        / "{assembly_id}"
        / "for-dramv"
        / "final-viral-combined-for-dramv.fa",
        tsv=VIRSORTER2
        / "{assembly_id}"
        / "for-dramv"
        / "viral-affi-contigs-for-dramv.tab",
        workdir=directory(VIRSORTER2 / "{assembly_id}"),
    log:
        VIRSORTER2 / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        workdir=lambda w: VIRSORTER2 / f"{w.assembly_id}",
    shell:
        """
        virsorter run \
            --working-dir {output.workdir} \
            --jobs {threads} \
            --prep-for-dramv \
            --tmpdir {output.workdir}/tmp \
            --rm-tmpdir \
            --verbose \
            --seqfile {input.fna} \
            --db-dir {input.database} \
        2> {log} 1>&2
        """


rule viral__virsorter:
    input:
        [
            VIRSORTER2
            / f"{assembly_id}"
            / "for-dramv"
            / "final-viral-combined-for-dramv.fa"
            for assembly_id in ASSEMBLIES
        ],
