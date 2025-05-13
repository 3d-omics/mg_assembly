rule prokaryotes__annotate__gtdbtk__classify_wf:
    """Run GTDB-Tk over the dereplicated genomes."""
    input:
        mags=PROK_MAGS,
        database=features["databases"]["gtdbtk"],
    output:
        work_dir=directory(PROK_GTDBTK),
    log:
        PROK_ANN / "gtdbtk.log",
    conda:
        "../../../environments/gtdbtk.yml"
    threads: 24
    resources:
        mem_mb=128 * 1024,
        runtime=24 * 60,
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


rule prokaryotes__annotate__gtdbtk__join_bac_and_ar:
    input:
        work_dir=PROK_GTDBTK,
    output:
        summary=PROK_ANN / "gtdbtk.summary.tsv",
        bac_tree=PROK_ANN / "gtdbtk.backbone.bac120.classify.tree",
        ar_tree=touch(PROK_ANN / "gtdbtk.ar53.classify.tree"),
    log:
        PROK_ANN / "gtdbtk.join.log",
    conda:
        "../../../environments/gtdbtk.yml"
    shell:
        """
        csvtk concat \
            --tabs \
            {input.work_dir}/gtdbtk.*.summary.tsv \
        > {output.summary} \
        2> {log}

        cp \
            --verbose \
            {input.work_dir}/classify/gtdbtk.backbone.bac120.classify.tree \
            {output.bac_tree} \
        2>> {log} 1>&2

        if [[ -f {input.work_dir}/classify/gtdbtk.ar53.classify.tree ]] ; then
            cp \
                --verbose \
                {input.work_dir}/classify/gtdbtk.ar53.classify.tree \
                {output.ar_tree} \
            2>> {log} 1>&2
        fi
        """


rule prokaryotes__annotate__gtdbtk__all:
    input:
        rules.prokaryotes__annotate__gtdbtk__join_bac_and_ar.output,
