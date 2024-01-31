rule _viral__cluster__genomad:
    input:
        fasta=MEGAHIT / "{assembly_id}.fa.gz",
        database=features["databases"]["genomad"],
    output:
        fna=GENOMADC / "{assembly_id}_summary" / "{assembly_id}_virus.fna",
        genes_tsv=GENOMADC / "{assembly_id}_summary" / "{assembly_id}_virus_genes.tsv",
        proteins=GENOMADC / "{assembly_id}_summary" / "{assembly_id}_virus_proteins.faa",
        summary_tsv=GENOMADC
        / "{assembly_id}_summary"
        / "{assembly_id}_virus_summary.tsv",
    log:
        GENOMADC / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        workdir=GENOMADC,
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


rule viral__cluster__genomad:
    input:
        [
            GENOMADC / f"{assembly_id}_summary" / f"{assembly_id}_virus.fna"
            for assembly_id in ASSEMBLIES
        ],
