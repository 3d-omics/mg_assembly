MVP_VOTU_FILTRATION = params["viral"]["mvp"]["votu_table"]["filtration"]
MVP_VOTU_NORMALIZATION = params["viral"]["mvp"]["votu_table"]["normalization"]


rule viruses__mvp__05_votu_table:
    """Merge all per-sample CoverM tables into vOTU tables (cohort-wide)"""
    input:
        metadata=VIR_MVP / "metadata.tsv",
        coverm=[
            VIR_MVP_READ_MAPPING / f"{sample_id}" / f"{sample_id}_CoverM.tsv"
            for sample_id in SAMPLES
        ],
    output:
        VIR_MVP_VOTU_TABLES
        / f"MVP_05_All_Sample_Filtered_{MVP_VOTU_FILTRATION}_Representative_Virus_Proviruses_vOTU_{MVP_VOTU_NORMALIZATION}_Table.tsv",
    log:
        VIR_MVP_VOTU_TABLES / "votu_table.log",
    conda:
        ENVS / "mvp.yml"
    params:
        covered_fraction=params["viral"]["mvp"]["votu_table"]["covered_fraction"],
        normalization=MVP_VOTU_NORMALIZATION,
        filtration=MVP_VOTU_FILTRATION,
        viral_min_genes=params["viral"]["mvp"]["votu_table"]["viral_min_genes"],
        host_viral_genes_ratio=params["viral"]["mvp"]["votu_table"][
            "host_viral_genes_ratio"
        ],
    shell:
        """
        mvip MVP_05_create_vOTU_table \
            --input {VIR_MVP} \
            --metadata {input.metadata} \
            --covered_fraction {params.covered_fraction} \
            --normalization {params.normalization} \
            --filtration {params.filtration} \
            --viral_min_genes {params.viral_min_genes} \
            --host_viral_genes_ratio {params.host_viral_genes_ratio} \
        2> {log} 1>&2
        """
