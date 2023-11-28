rule _dereplicate__gtdbtk__download:
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


rule _dereplicate__gtdbtk__classify:
    """Run GTDB-Tk over the dereplicated genomes."""
    input:
        bin_folder=DREP / "dereplicated_genomes",
        database=features["gtdbtk_database"],
    output:
        summary=DREP_GTDBTK / "gtdbtk.summary.tsv",
        align=DREP_GTDBTK / "align.tar.gz",
        classify=DREP_GTDBTK / "classify.tar.gz",
        identify=DREP_GTDBTK / "identify.tar.gz",
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
        mem_mb=64 * 1024,
        runtime=24 * 60,
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
            cat {params.bac120} > {output.summary} 2>> {log}
        fi

        for folder in align classify identify ; do
            tar \
                --create \
                --directory {params.out_dir} \
                --file {params.out_dir}/${{folder}}.tar.gz \
                --remove-files \
                --use-compress-program="pigz --processes {threads}" \
                --verbose \
                ${{folder}} \
            2>> {log} 1>&2
        done

        rm \
            --recursive \
            --force \
            --verbose \
            {params.ar53} \
            {params.bac120} \
        2>> {log} 1>&2
        """


rule dereplicate__gtdbtk:
    """Run the gtdbtk subworkflow"""
    input:
        rules._dereplicate__gtdbtk__classify.output,
