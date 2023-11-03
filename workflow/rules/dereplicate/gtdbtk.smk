rule dereplicate_gtdbtk_classify:
    """Run GTDB-Tk over the dereplicated genomes.""" ""
    input:
        bin_folder=DREP / "dereplicated_genomes",
        database=features["gtdbtk_database"],
    output:
        ar53=touch(DREP_GTDBTK / "gtdbtk.ar53.summary.tsv"),
        bac120=touch(DREP_GTDBTK / "gtdbtk.bac120.summary.tsv"),
    log:
        DREP_GTDBTK / "gtdbtk_classify.log",
    conda:
        "gtdbtk.yml"
    params:
        out_dir=DREP_GTDBTK,
    threads: 24
    resources:
        mem_mb=double_ram(params["dereplicate"]["gtdbtk"]["memory_gb"]),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        export GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {input.bin_folder} \
            --extension fa \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --skip_ani_screen \
        2> {log} 1>&2
        """


rule dereplicate_gtdbtk_join_summaries:
    """Join both bacterial and archaeal summaries"""
    input:
        ar53=DREP_GTDBTK / "gtdbtk.ar53.summary.tsv",
        bac120=DREP_GTDBTK / "gtdbtk.bac120.summary.tsv",
    output:
        summary=DREP_GTDBTK / "gtdbtk.summary.tsv",
    log:
        DREP_GTDBTK / "gtdbtk.summary.log",
    conda:
        "gtdbtk.yml"
    shell:
        """
        ( csvstack \
            --tabs \
            {input.bac120} \
            {input.ar53} \
        | csvformat \
            --out-tabs \
        > {output.summary} ) \
        2> {log}
        """


rule dereplicate_gtdbtk:
    """Run the gtdbtk subworkflow"""
    input:
        rules.dereplicate_gtdbtk_join_summaries.output,
