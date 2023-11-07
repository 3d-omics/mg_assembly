rule dereplicate_checkm2_download:
    output:
        features["checkm2_database"],
    log:
        DREP_CHECKM / "download.log",
    conda:
        "checkm2.yml"
    shell:
        """
        checkm2 database \
            --download \
            --path $(dirname {output}) \
        2> {log} 1>&2

        mv $(dirname {output})/CheckM2_database/* {output} 1>&2
        rmdir $(dirname {output})/CheckM2_database 2>>{log} 1>&2
        """


rule dereplicate_checkm2_predict:
    """Run CheckM2 over the dereplicated mags"""
    input:
        mags=DREP / "dereplicated_genomes",
        db=features["checkm2_database"],
    output:
        DREP_CHECKM / "quality_report.tsv",
    log:
        DREP_CHECKM / "quality_report.log",
    threads: 24
    conda:
        "checkm2.yml"
    params:
        out_dir=DREP_CHECKM / "predict",
    shell:
        """
        checkm2 predict \
            --threads {threads} \
            --input {input.mags} \
            --extension .fa \
            --output-directory {params.out_dir} \
            --database_path {input.db} \
            --remove_intermediates \
        2> {log} 1>&2

        mv {params.out_dir}/quality_report.tsv {output} 2>> {log} 1>&2
        rm --recursive --verbose --force {params.out_dir} 2>> {log} 1>&2
        """


rule dereplicate_checkm2:
    """Run CheckM2"""
    input:
        rules.dereplicate_checkm2_predict.output,
