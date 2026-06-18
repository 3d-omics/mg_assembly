# NOTE: vRhyme writes one shared "07_BINNING/07A_vRHYME_OUTPUT" directory regardless
# of --binning_sample_group: concurrent per-sample invocations would race on the same
# directory, so unlike 01/02/04 this stage cannot be parallelized per sample.


rule viruses__mvp__07_binning:
    """vRhyme binning of the representative viral catalog using all samples' read mapping"""
    input:
        metadata=VIR_MVP / "metadata.tsv",
        representative_fasta=VIR_MVP_CLUSTERING
        / "MVP_03_All_Sample_Filtered_Relaxed_Representative_Virus_Provirus_Sequences.fna",
        bams=[
            VIR_MVP_READ_MAPPING / f"{sample_id}" / f"{sample_id}_sorted.bam"
            for sample_id in SAMPLES
        ],
    output:
        directory(VIR_MVP_BINNING),
    log:
        VIR_MVP / "mvp_07_binning.log",
    conda:
        ENVS / "mvp.yml"
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=24 * 60,
    params:
        normalization=params["viral"]["mvp"]["binning"]["normalization"],
        filtration=params["viral"]["mvp"]["binning"]["filtration"],
        viral_min_genes=params["viral"]["mvp"]["binning"]["viral_min_genes"],
        host_viral_genes_ratio=params["viral"]["mvp"]["binning"][
            "host_viral_genes_ratio"
        ],
        read_type=params["viral"]["mvp"]["binning"]["read_type"],
        interleaved=params["viral"]["mvp"]["binning"]["interleaved"],
    shell:
        """
        mvip MVP_07_do_binning \
            --input {VIR_MVP} \
            --metadata {input.metadata} \
            --normalization {params.normalization} \
            --filtration {params.filtration} \
            --viral_min_genes {params.viral_min_genes} \
            --host_viral_genes_ratio {params.host_viral_genes_ratio} \
            --read_type {params.read_type} \
            --interleaved {params.interleaved} \
            --threads {threads} \
        2> {log} 1>&2
        """
