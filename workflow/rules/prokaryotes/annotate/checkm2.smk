rule prokaryotes__annotate__checkm2__:
    """Run CheckM2 over the dereplicated mags"""
    input:
        mags=MAGS,
        db=features["databases"]["checkm2"],
    output:
        report=PROK_ANN / "checkm2.quality_report.tsv",
        tmp_dir = temp(directory(PROK_ANN / "checkm2.quality_report"))
    log:
        PROK_ANN / "quality_report.log",
    conda:
        "checkm2.yml"
    shell:
        """
        checkm2 predict \
            --threads {threads} \
            --input {input.mags} \
            --extension .fa \
            --output-directory {output.tmp_dir} \
            --database_path {input.db}/uniref100.KO.1.dmnd \
            --remove_intermediates \
        2>> {log} 1>&2

        cp \
            --verbose \
            {output.tmp_dir}/quality_report.tsv \
            {output.report} \
        2>> {log} 1>&2
        """


rule prokaryotes__annotate__checkm2:
    input:
        rules.prokaryotes__annotate__checkm2__.output,
