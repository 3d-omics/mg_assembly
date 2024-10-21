rule viruses__cluster__genomad__:
    input:
        fasta=ASSEMBLE_MEGAHIT / "{assembly_id}.fa.gz",
        database=features["databases"]["genomad"],
    output:
        plasmid=GENOMADC / "{assembly_id}_plasmid.fna.gz",
        plasmid_genes=GENOMADC / "{assembly_id}_plasmid_genes.tsv.gz",
        plasmid_proteins=GENOMADC / "{assembly_id}_plasmid_proteins.faa.gz",
        plasmid_summary=GENOMADC / "{assembly_id}_plasmid_summary.tsv.gz",
        json=GENOMADC / "{assembly_id}_summary.json.gz",
        virus=GENOMADC / "{assembly_id}_virus.fna.gz",
        virus_genes=GENOMADC / "{assembly_id}_virus_genes.tsv.gz",
        virus_proteins=GENOMADC / "{assembly_id}_virus_proteins.faa.gz",
        virus_summary_tsv=GENOMADC / "{assembly_id}_virus_summary.tsv.gz",
    log:
        GENOMADC / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        genomad_workdir=GENOMADC,
        genomad_summary_dir=lambda w: GENOMADC / f"{w.assembly_id}_summary",
        extra=params["viral"]["genomad"]["extra"],
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
            {params.genomad_workdir} \
            {input.database} \
        2> {log} 1>&2

        bgzip \
            --threads {threads} \
            {params.genomad_summary_dir}/* \
        2>> {log}

        mv \
            --verbose \
            {params.genomad_summary_dir}/* \
            {GENOMADC} \
        2>> {log} 1>&2
        """


rule viruses__cluster__genomad:
    input:
        [GENOMADC / f"{assembly_id}_virus.fna.gz" for assembly_id in ASSEMBLIES],
