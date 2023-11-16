rule dereplicate_gtdbtk_download:
    output:
        directory(features["gtdbtk_database"]),
    log:
        features["gtdbtk_database"] + ".log",
    conda:
        "gtdbtk.yml"
    shell:
        """
        mkdir --parents {output} 2> {log} 1>&2

        wget \
            --continue \
            --directory {output} \
            https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz \
        2>> {log} 1>&2

        tar \
            --extract \
            --verbose \
            --file {output}/gtdbtk_data.tar.gz \
            --directory {output} \
        2>> {log} 1>&2

        mv \
            {output}/release*/* \
            {output}/ \
        2>> {log} 1>&2
        """


rule dereplicate_gtdbtk_classify:
    """Run GTDB-Tk over the dereplicated genomes."""
    input:
        bin_folder=DREP / "dereplicated_genomes",
        database=features["gtdbtk_database"],
    output:
        summary=DREP_GTDBTK / "gtdbtk.summary.tsv",
    log:
        DREP_GTDBTK / "gtdbtk_classify.log",
    conda:
        "gtdbtk.yml"
    params:
        out_dir=DREP_GTDBTK,
        ar53=DREP_GTDBTK / "gtdbtk.ar53.summary.tsv",
        bac120=DREP_GTDBTK / "gtdbtk.bac120.summary.tsv",
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

        if [[ -f {params.ar53} ]] ; then
            ( csvstack \
                --tabs \
                {params.bac120} \
                {params.ar53} \
            | csvformat \
                --out-tabs \
            > {output.summary} ) \
            2>> {log}
        else
            cp {params.bac120} {output.summary} 2>> {log}
        fi

        for folder in align classify identify ; do
            tar \
                --create \
                --verbose \
                --remove-files \
                --use-compress-program="pigz --processes {threads}" \
                --file {params.out_dir}/${{folder}}.tar.gz \
                {params.out_dir}/${{folder}} \
            2>> {log} 1>&2
        done
        """


rule dereplicate_gtdbtk:
    """Run the gtdbtk subworkflow"""
    input:
        rules.dereplicate_gtdbtk_classify.output,
