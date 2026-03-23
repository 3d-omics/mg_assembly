rule viruses__cluster__genomad:
    input:
        fasta=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
        database=features["databases"]["genomad"],
    output:
        plasmid=VIR_GENOMADC / "{assembly_id}_plasmid.fna.gz",
        plasmid_genes=VIR_GENOMADC / "{assembly_id}_plasmid_genes.tsv.gz",
        plasmid_proteins=VIR_GENOMADC / "{assembly_id}_plasmid_proteins.faa.gz",
        plasmid_summary=VIR_GENOMADC / "{assembly_id}_plasmid_summary.tsv.gz",
        json=VIR_GENOMADC / "{assembly_id}_summary.json.gz",
        virus=VIR_GENOMADC / "{assembly_id}_virus.fna.gz",
        virus_genes=VIR_GENOMADC / "{assembly_id}_virus_genes.tsv.gz",
        virus_proteins=VIR_GENOMADC / "{assembly_id}_virus_proteins.faa.gz",
        virus_summary_tsv=VIR_GENOMADC / "{assembly_id}_virus_summary.tsv.gz",
    log:
        VIR_GENOMADC / "{assembly_id}.log",
    conda:
        "../../../environments/genomad.yml"
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        genomad_workdir=VIR_GENOMADC,
        genomad_summary_dir=lambda w: VIR_GENOMADC / f"{w.assembly_id}_summary",
        extra=params["viral"]["genomad"]["extra"],
        use_cuda=params["viral"]["genomad"]["use_cuda"],
    shadow:
        "minimal"
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
    shell:
        """
        if [[ $(gzip -dc {input.fasta} | wc -l ) -lt 2 ]] ; then
            echo "Empty fasta. Touching outputs" 2> {log}  1>&2
            touch {output}
            exit 0
        fi

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
            {params.genomad_workdir} \
        2>> {log} 1>&2
        """


rule viruses__cluster__genomad__all:
    input:
        [VIR_GENOMADC / f"{assembly_id}_virus.fna.gz" for assembly_id in ASSEMBLIES],
