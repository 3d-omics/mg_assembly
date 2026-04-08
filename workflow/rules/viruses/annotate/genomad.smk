rule viruses__annotate__genomad:
    input:
        fasta=VIR_CLUSTER / "mmseqs.rep_seq.fa.gz",
        database=features["databases"]["genomad"],
    output:
        plasmid=VIR_GENOMADA / "rep_seq_plasmid.fna.gz",
        plasmid_genes=VIR_GENOMADA / "rep_seq_plasmid_genes.tsv.gz",
        plasmid_proteins=VIR_GENOMADA / "rep_seq_plasmid_proteins.faa.gz",
        plasmid_summary=VIR_GENOMADA / "rep_seq_plasmid_summary.tsv.gz",
        json=VIR_GENOMADA / "rep_seq_summary.json.gz",
        virus=VIR_GENOMADA / "rep_seq_virus.fna.gz",
        virus_genes=VIR_GENOMADA / "rep_seq_virus_genes.tsv.gz",
        virus_proteins=VIR_GENOMADA / "rep_seq_virus_proteins.faa.gz",
        virus_summary=VIR_GENOMADA / "rep_seq_virus_summary.tsv.gz",
    log:
        VIR_GENOMADA / "genomad.log",
    conda:
        "../../../environments/genomad.yml"
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        workdir=VIR_GENOMADA,
        extra=params["viral"]["genomad"]["extra"],
        tmp_prefix=VIR_GENOMADA / "rep_seq_summary",
        use_cuda=params["viral"]["genomad"]["use_cuda"],
    shadow:
        "minimal"
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=60,
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


rule viruses__annotate__genomad__all:
    input:
        rules.viruses__annotate__genomad.output,
