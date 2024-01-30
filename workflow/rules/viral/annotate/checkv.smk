rule _viral__annotate__checkv:
    input:
        fna=MMSEQS / "results_all_seqs.fasta",
        database=features["databases"]["checkv"],
    output:
        fna=CHECKVA / "all.fna",
        summary=CHECKVA / "quality_summary.tsv",
    log:
        CHECKVA / "checkv.log",
    conda:
        "__environment__.yml"
    params:
        workdir=CHECKVA,
    threads: 8
    shell:
        """
        checkv end_to_end \
            {input.fna} \
            {params.workdir} \
            --remove_tmp \
            -t {threads} \
            --restart \
            -d {input.database} \
        2> {log} 1>&2

        cat \
            {params.workdir}/proviruses.fna \
            {params.workdir}/viruses.fna \
        > {output.fna} \
        2>> {log}
        """


rule viral__annotate__checkv:
    input:
        rules._viral__annotate__checkv.input,
