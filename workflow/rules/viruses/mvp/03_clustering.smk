rule viruses__mvp__03_clustering:
    """Cohort-wide ANI clustering of filtered viral sequences; also builds the MVP_04 bowtie2 index"""
    input:
        metadata=VIR_MVP / "metadata.tsv",
        fastas=[
            VIR_MVP_CHECKV
            / f"{sample_id}"
            / f"MVP_02_{sample_id}_Filtered_Relaxed_Virus_Provirus_Sequences.fna"
            for sample_id in SAMPLES
        ],
    output:
        representative_fasta=VIR_MVP_CLUSTERING
        / "MVP_03_All_Sample_Filtered_Relaxed_Representative_Virus_Provirus_Sequences.fna",
        representative_summary=VIR_MVP_CLUSTERING
        / "MVP_03_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Quality_Summary.tsv",
        clusters=VIR_MVP_CLUSTERING
        / "MVP_03_All_Sample_Filtered_Relaxed_Virus_Provirus_Sequences_Clustering_ANI_Clusters.tsv",
        bowtie2_index=multiext(
            str(VIR_MVP_READ_MAPPING / "reference"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        VIR_MVP_CLUSTERING / "clustering.log",
    conda:
        ENVS / "mvp.yml"
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=12 * 60,
    params:
        min_ani=params["viral"]["mvp"]["clustering"]["min_ani"],
        min_tcov=params["viral"]["mvp"]["clustering"]["min_tcov"],
        min_qcov=params["viral"]["mvp"]["clustering"]["min_qcov"],
        read_type=params["viral"]["mvp"]["clustering"]["read_type"],
    shell:
        """
        mvip MVP_03_do_clustering \
            --input {VIR_MVP} \
            --metadata {input.metadata} \
            --min_ani {params.min_ani} \
            --min_tcov {params.min_tcov} \
            --min_qcov {params.min_qcov} \
            --read-type {params.read_type} \
            --threads {threads} \
        2> {log} 1>&2
        """
