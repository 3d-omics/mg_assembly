rule dereplicate_gtdbtk:
    """Run GTDB-Tk over the dereplicated genomes.""" ""
    input:
        bin_folder=DREP / "dereplicated_genomes",
        database=features["gtdbtk_database"],
    output:
        outdir=directory(DREP_GTDBTK),
    log:
        DREP_GTDBTK / "gtdbtk.log",
    conda:
        "gtdbtk.yml"
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
            --out_dir {output.outdir} \
            --cpus {threads} \
            --skip_ani_screen \
            --full_tree \
        2> {log} 1>&2
        """
