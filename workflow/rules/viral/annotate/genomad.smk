rule _viral__annotate__genomad:
    input:
        fasta=MMSEQS / "rep_seq.fasta",
        database=features["databases"]["genomad"],
    output:
        plasmid=GENOMADA  / "rep_seq_plasmid.fna",
        plasmid_genes=GENOMADA / "rep_seq_plasmid_genes.tsv",
        plasmid_proteins=GENOMADA / "rep_seq_plasmid_proteins.faa",
        plasmid_summary=GENOMADA / "rep_seq_plasmid_summary.tsv",
        json=GENOMADA / "rep_seq_summary.json",
        virus=GENOMADA / "rep_seq_virus.fna",
        virus_genes=GENOMADA / "rep_seq_virus_genes.tsv",
        virus_proteins=GENOMADA / "rep_seq_virus_proteins.faa",
        virus_summary=GENOMADA / "rep_seq_virus_summary.tsv",
    log:
        GENOMADA / "genomad.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        workdir=GENOMADA,
        extra=params["viral"]["genomad"]["extra"],
        tmp_prefix=GENOMADA / "rep_seq_summary",
        use_cuda=params["viral"]["genomad"]["use_cuda"],
    resources:
        mem_mb=32 * 1024,
    shadow: "minimal",
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
