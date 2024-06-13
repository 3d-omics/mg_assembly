rule prokaryotes__annotate__gtdbtk__:
    """Run GTDB-Tk over the dereplicated genomes."""
    input:
        fasta_folder=DREP / "dereplicated_genomes",
        database=features["databases"]["gtdbtk"],
    output:
        summary=GTDBTK / "gtdbtk.summary.tsv",
        align=GTDBTK / "align.tar.gz",
        classify=GTDBTK / "classify.tar.gz",
        identify=GTDBTK / "identify.tar.gz",
        bac_tree=GTDBTK / "gtdbtk.backbone.bac120.classify.tree",
        ar_tree=touch(GTDBTK / "gtdbtk.ar53.tree"),
    log:
        GTDBTK / "gtdbtk_classify.log",
    conda:
        "gtdbtk.yml"
    params:
        out_dir=GTDBTK,
        ar53=GTDBTK / "gtdbtk.ar53.summary.tsv",
        bac120=GTDBTK / "gtdbtk.bac120.summary.tsv",
        ar_tree=GTDBTK / "classify" / "",
        bac_tree=GTDBTK / "classify" / "gtdbtk.backbone.bac120.classify.tree",
    resources:
        attempt=get_attempt,
    retries: 5
    shell:
        """
        rm \
            --recursive \
            --force \
            {params.out_dir}/align \
            {params.out_dir}/classify \
            {params.out_dir}/identify \
            {params.out_dir}/gtdbtk_classify.log \
            {params.out_dir}/gtdbtk.json \
            {params.out_dir}/gtdbtk.summary.tsv \
            {params.out_dir}/gtdbtk.warnings.log \
            {params.ar53} \
            {params.bac120} \
        2>> {log}.{resources.attempt} 1>&2

        export GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {input.fasta_folder} \
            --extension fa.gz \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --skip_ani_screen \
        2>> {log}.{resources.attempt} 1>&2

        if [[ -f {params.ar53} ]] ; then
            ( csvstack \
                --tabs \
                {params.bac120} \
                {params.ar53} \
            | csvformat \
                --out-tabs \
            > {output.summary} \
            ) 2>> {log}.{resources.attempt}
        else
            cp {params.bac120} {output.summary} 2>> {log}.{resources.attempt}
        fi

        cp --verbose {params.bac_tree} {output.bac_tree} 2>> {log}

        if [[ -f {params.ar_tree} ]] ; then
            cp --verbose {params.ar_tree} {output.ar_tree} 2>> {log} 1>&2
        else
            touch {output.ar_tree} 2>> {log} 1>&2
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
            2>> {log}.{resources.attempt} 1>&2
        done

        mv {log}.{resources.attempt} {log}
        """


rule prokaryotes__annotate__gtdbtk:
    input:
        rules.prokaryotes__annotate__gtdbtk__.output,
