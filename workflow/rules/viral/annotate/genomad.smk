rule _viral__annotate__genomad:
    input:
        fasta=MMSEQS / "results_all_seqs.fasta",
        database=features["databases"]["genomad"],
    output:
        fna=GENOMADA / "all_summary" / "all_virus.fna",
        genes_tsv=GENOMADA / "all_summary" / "all_virus_genes.tsv",
        proteins=GENOMADA / "all_summary" / "all_virus_proteins.faa",
        summary_tsv=GENOMADA / "all_summary" / "all_virus_summary.tsv",
    log:
        GENOMADA / "all.log",
    conda:
        "__environment__.yml"
    threads: 8
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        workdir=GENOMADA,
        extra=params["viral"]["genomad"]["extra"],
    resources:
        mem_mb=32 * 1024,
    shell:
        """
        genomad end-to-end \
            {params.filtering} \
            --cleanup \
            --restart \
            --verbose \
            --threads {threads} \
            --disable-nn-classification \
            {params.extra} \
            {input.fasta} \
            {params.workdir} \
            {input.database} \
        2> {log} 1>&2
        """


rule viral__annotate__genomad:
    input:
        rules._viral__annotate__genomad.output,
