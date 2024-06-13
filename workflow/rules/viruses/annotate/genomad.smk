rule viruses__annotate__genomad__:
    input:
        fasta=MMSEQS / "rep_seq.fasta.gz",
        database=features["databases"]["genomad"],
    output:
        plasmid=GENOMADA / "rep_seq_plasmid.fna.gz",
        plasmid_genes=GENOMADA / "rep_seq_plasmid_genes.tsv.gz",
        plasmid_proteins=GENOMADA / "rep_seq_plasmid_proteins.faa.gz",
        plasmid_summary=GENOMADA / "rep_seq_plasmid_summary.tsv.gz",
        json=GENOMADA / "rep_seq_summary.json.gz",
        virus=GENOMADA / "rep_seq_virus.fna.gz",
        virus_genes=GENOMADA / "rep_seq_virus_genes.tsv.gz",
        virus_proteins=GENOMADA / "rep_seq_virus_proteins.faa.gz",
        virus_summary=GENOMADA / "rep_seq_virus_summary.tsv.gz",
    log:
        GENOMADA / "genomad.log",
    conda:
        "__environment__.yml"
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        workdir=GENOMADA,
        extra=params["viral"]["genomad"]["extra"],
        tmp_prefix=GENOMADA / "rep_seq_summary",
        use_cuda=params["viral"]["genomad"]["use_cuda"],
    shadow:
        "minimal"
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

        bgzip \
            --threads {threads} \
            {params.tmp_prefix}/* \
        2>> {log} 1>&2

        mv \
            {params.tmp_prefix}/* \
            {params.workdir} \
        2>> {log} 1>&2
        """


rule viruses__annotate__genomad:
    input:
        rules.viruses__annotate__genomad__.output,
