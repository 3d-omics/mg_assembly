rule _viral__annotate__genomad:
    input:
        fasta=MMSEQS / "cluster.fa",
        database=features["databases"]["genomad"],
    output:
        fna=GENOMADA / "cluster_virus.fna",
        genes_tsv=GENOMADA / "cluster_virus_genes.tsv",
        proteins=GENOMADA / "cluster_virus_proteins.faa",
        summary_tsv=GENOMADA / "cluster_virus_summary.tsv",
    log:
        GENOMADA / "genomad.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        workdir=GENOMADA,
        extra=params["viral"]["genomad"]["extra"],
        tmp_prefix=GENOMADA / "cluster_summary",
        use_cuda=params["viral"]["genomad"]["use_cuda"],
    resources:
        mem_mb=32 * 1024,
    shell:
        """
        {params.use_cuda}

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

        mv \
            {params.tmp_prefix}/* \
            {params.workdir} \
        2>> {log} 1>&2
        """


rule viral__annotate__genomad:
    input:
        rules._viral__annotate__genomad.output,
