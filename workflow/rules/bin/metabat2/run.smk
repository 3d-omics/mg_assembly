include: "functions.smk"


rule metabat2_prepare_one:
    input:
        bams=get_bams_for_metabat2,
    output:
        depth=METABAT2 / "prepare" / "{assembly_id}.depth",
        paired_contigs=METABAT2 / "prepare" / "{assembly_id}.paired",
    log:
        METABAT2 / "prepare" / "{assembly_id}.log",
    conda:
        "metabat2.yml"
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            --pairedContigs {output.paired_contigs} \
            {input.bams} \
        2> {log} 1>&2
        """


rule metabat2_run_one:
    input:
        assembly=ASSEMBLE_RENAME / "{assembly_id}.fa",
        depth=METABAT2 / "prepare/{assembly_id}.depth",
    output:
        bins=directory(METABAT2 / "bins/{assembly_id}/"),
    log:
        METABAT2 / "bins" / "{assembly_id}.log",
    conda:
        "metabat2.yml"
    params:
        bins_prefix=lambda wildcards: METABAT2
        / f"bins/{wildcards.assembly_id}/{wildcards.assembly_id}",
    threads: 24
    shell:
        """
        metabat2 \
            --inFile {input.assembly} \
            --abdFile {input.depth} \
            --outFile {params.bins_prefix} \
            --numThreads {threads} \
        2> {log} 1>&2
        """


rule metabat2:
    input:
        [METABAT2 / f"bins/{assembly_id}" for assembly_id in ASSEMBLIES],
