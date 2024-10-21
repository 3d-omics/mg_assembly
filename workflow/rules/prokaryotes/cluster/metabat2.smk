rule prokaryotes__cluster__metabat2__:
    """Run metabat2 end-to-end on a single assembly"""
    input:
        bams=get_bams_from_assembly_id,
        assembly=ASSEMBLE_MEGAHIT / "{assembly_id}.fa.gz",
    output:
        bins=directory(METABAT2 / "{assembly_id}"),
    log:
        METABAT2 / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    params:
        bins_prefix=lambda w: METABAT2 / f"{w.assembly_id}/bin",
        depth=lambda w: METABAT2 / f"{w.assembly_id}.depth",
        paired=lambda w: METABAT2 / f"{w.assembly_id}.paired",
        workdir=METABAT2,
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {params.depth} \
            --pairedContigs {params.paired} \
            {input.bams} \
        2>> {log} 1>&2

        metabat2 \
            --inFile {input.assembly} \
            --abdFile {params.depth} \
            --outFile {params.bins_prefix} \
            --numThreads {threads} \
            --verbose \
        2> {log} 1>&2

        rm \
            --force \
            --verbose \
            {params.depth} \
            {params.paired} \
        2>> {log} 1>&2

        parallel --jobs {threads} \
            bgzip \
                --compress-level 9 \
                {{}} \
        ::: {output.bins}/*.fa \
        2>> {log} 1>&2
        """


rule prokaryotes__cluster__metabat2:
    """Run metabat2 over all assemblies"""
    input:
        [METABAT2 / assembly_id for assembly_id in ASSEMBLIES],
