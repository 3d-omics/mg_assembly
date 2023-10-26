rule bin_metabat2_prepare_one:
    """Compute coverages for metabat2"""
    input:
        bams=get_bams_from_assembly_id,
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


rule bin_metabat2_run_one:
    """Run metabat2 over a single assembly"""
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
        bins_prefix=compose_bins_prefix_for_metabat2_run_one,
    threads: 24
    resources:
        runtime=24 * 60,
        mem_mb=8 * 1024,
    shell:
        """
        metabat2 \
            --inFile {input.assembly} \
            --abdFile {input.depth} \
            --outFile {params.bins_prefix} \
            --numThreads {threads} \
        2> {log} 1>&2
        """


rule bin_metabat2:
    """Run metabat2 over all assemblies"""
    input:
        [METABAT2 / f"bins/{assembly_id}" for assembly_id in ASSEMBLIES],
