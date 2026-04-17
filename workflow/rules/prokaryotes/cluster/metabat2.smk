rule prokaryotes__cluster__metabat2:
    """Run metabat2 end-to-end on a single assembly"""
    input:
        bams=get_bams_from_assembly_id,
        assembly=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
    output:
        bins=directory(PROK_METABAT2 / "{assembly_id}"),
    log:
        PROK_METABAT2 / "{assembly_id}.log",
    conda:
        ENVS / "metabat2.yml"
    threads: 24
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=24 * 60,
    params:
        bins_prefix=lambda w: PROK_METABAT2 / w.assembly_id / "bin",
        depth=lambda w: PROK_METABAT2 / f"{w.assembly_id}.depth",
        paired=lambda w: PROK_METABAT2 / f"{w.assembly_id}.paired",
        workdir=PROK_METABAT2,
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
                {{}} \
        ::: {output.bins}/*.fa \
        2>> {log} 1>&2
        """


rule prokaryotes__cluster__metabat2__all:
    """Run metabat2 over all assemblies"""
    input:
        [PROK_METABAT2 / assembly_id for assembly_id in ASSEMBLIES],
