rule viruses__mvp__02_filter:
    """Merge and filter geNomad+CheckV outputs for one sample, in parallel across samples"""
    input:
        metadata=VIR_MVP / "metadata.tsv",
        summary=VIR_MVP_CHECKV
        / "{sample_id}"
        / "MVP_01_{sample_id}_Unfiltered_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv",
    output:
        summary=VIR_MVP_CHECKV
        / "{sample_id}"
        / "MVP_02_{sample_id}_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv",
        fasta=VIR_MVP_CHECKV
        / "{sample_id}"
        / "MVP_02_{sample_id}_Filtered_Relaxed_Virus_Provirus_Sequences.fna",
    log:
        VIR_MVP_CHECKV / "{sample_id}.filter.log",
    conda:
        ENVS / "mvp.yml"
    params:
        sample_number=get_sample_number,
        viral_min_genes=params["viral"]["mvp"]["filter"]["viral_min_genes"],
        host_viral_genes_ratio=params["viral"]["mvp"]["filter"][
            "host_viral_genes_ratio"
        ],
    shell:
        """
        mvip MVP_02_filter_genomad_checkv \
            --input {VIR_MVP} \
            --metadata {input.metadata} \
            --sample_group {params.sample_number} \
            --viral_min_genes {params.viral_min_genes} \
            --host_viral_genes_ratio {params.host_viral_genes_ratio} \
        2> {log} 1>&2
        """


rule viruses__mvp__02_filter__all:
    input:
        [
            VIR_MVP_CHECKV
            / f"{sample_id}"
            / f"MVP_02_{sample_id}_Filtered_Relaxed_Virus_Provirus_Sequences.fna"
            for sample_id in SAMPLES
        ],
