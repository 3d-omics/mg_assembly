rule viruses__mvp__01_genomad_checkv:
    """Run geNomad + CheckV for one sample (identification), in parallel across samples"""
    input:
        metadata=VIR_MVP / "metadata.tsv",
        setup=VIR_MVP / "MVP_00_Summary_Report.txt",
        genomad_db=features["databases"]["genomad"],
        checkv_db=features["databases"]["checkv"],
    output:
        summary=VIR_MVP_CHECKV
        / "{sample_id}"
        / "MVP_01_{sample_id}_Unfiltered_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv",
    log:
        VIR_MVP_GENOMAD / "{sample_id}.log",
    conda:
        ENVS / "mvp.yml"
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=12 * 60,
    params:
        sample_number=get_sample_number,
        genomad_mode=params["viral"]["mvp"]["genomad_checkv"]["genomad_mode"],
        min_seq_size=params["viral"]["mvp"]["genomad_checkv"]["min_seq_size"],
    shell:
        """
        mvip MVP_01_run_genomad_checkv \
            --input {VIR_MVP} \
            --metadata {input.metadata} \
            --sample_group {params.sample_number} \
            --min_seq_size {params.min_seq_size} \
            {params.genomad_mode} \
            --genomad_db_path {input.genomad_db} \
            --checkv_db_path {input.checkv_db} \
            --threads {threads} \
        2> {log} 1>&2
        """


rule viruses__mvp__01_genomad_checkv__all:
    input:
        [
            VIR_MVP_CHECKV
            / f"{sample_id}"
            / f"MVP_01_{sample_id}_Unfiltered_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv"
            for sample_id in SAMPLES
        ],
