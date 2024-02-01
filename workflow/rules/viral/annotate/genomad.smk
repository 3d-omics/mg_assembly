rule _viral__annotate__genomad:
    input:
        fasta=MMSEQS / "results_all_seqs.fasta",
        database=features["databases"]["genomad"],
    output:
        fna=GENOMADA / "results_all_seqs_summary" / "all_virus.fna",
        genes_tsv=GENOMADA / "results_all_seqs_all_summary" / "results_all_seqs_virus_genes.tsv",
        proteins=GENOMADA / "results_all_seqs_all_summary" / "results_all_seqs_virus_proteins.faa",
        summary_tsv=GENOMADA / "results_all_seqs_all_summary" / "results_all_seqs_virus_summary.tsv",
    log:
        GENOMADA / "all.log",
    conda:
        "__environment__.yml"
    threads: 24
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
            {params.extra} \
            {input.fasta} \
            {params.workdir} \
            {input.database} \
        2> {log} 1>&2
        """


rule viral__annotate__genomad:
    input:
        rules._viral__annotate__genomad.output,
