rule prokaryotes__annotate__gtdbtk__classify_wf__:
    """Run GTDB-Tk over the dereplicated genomes."""
    input:
        mags=MAGS,
        database=features["databases"]["gtdbtk"],
    output:
        work_dir=directory(GTDBTK),
    log:
        PROK_ANN / "gtdbtk.log",
    conda:
        "__environment__.yml"
    shell:
        """
        export GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {input.mags} \
            --extension fa \
            --out_dir {output.work_dir} \
            --cpus {threads} \
            --skip_ani_screen \
        2> {log} 1>&2
        """


rule prokaryotes__annotate__gtdbtk__join_bac_and_ar__:
    input:
        work_dir=GTDBTK,
    output:
        summary=PROK_ANN / "gtdbtk.summary.tsv",
        bac_tree=PROK_ANN / "gtdbtk.backbone.bac120.classify.tree",
        ar_tree=touch(PROK_ANN / "gtdbtk.ar53.classify.tree"),
    log:
        PROK_ANN / "gtdbtk.join.log",
    conda:
        "__environment__.yml"
    shell:
        """
        if [[ -f {input.work_dir}/gtdbtk.ar122.summary.tsv ]] ; then

            csvstack \
                --tabs \
                {input.work_dir}/gtdbtk.bac120.summary.tsv \
                {input.work_dir}/gtdbtk.ar53.summary.tsv \
            | csvformat \
                --out-tabs \
            > {output.summary} \
            2> {log}

            cp \
                --verbose \
                {input.work_dir}/classify/gtdbtk.ar53.classify.tree \
                {output.ar_tree} \
            2>> {log}

        else

            cp \
                --verbose \
                {input.work_dir}/gtdbtk.bac120.summary.tsv \
                {output.summary} \

        fi

        cp \
            --verbose \
            {input.work_dir}/classify/gtdbtk.backbone.bac120.classify.tree \
            {output.bac_tree} \
        2>> {log} 1>&2
        """


rule prokaryotes__annotate__gtdbtk:
    input:
        rules.prokaryotes__annotate__gtdbtk__join_bac_and_ar__.output,
