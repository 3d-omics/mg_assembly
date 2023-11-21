rule _dereplicate__checkm2__download:
    output:
        directory(features["checkm2_database"]),
    log:
        DREP_CHECKM / "download.log",
    conda:
        "checkm2.yml"
    shell:
        """
        mkdir --parents {output} 2> {log} 1>&2

        checkm2 database \
            --download \
            --path $(dirname {output}) \
        2>> {log} 1>&2

        mv $(dirname {output})/CheckM2_database/* {output}/ 2>> {log} 1>&2
        rmdir $(dirname {output})/CheckM2_database 2>>{log} 1>&2
        """


rule _dereplicate__checkm2__predict:
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
    resources:
        mem_mb=16 * 1024,
    shell:
        """
        rm -rfv {params.out_dir} 2> {log} 1>&2

        checkm2 predict \
            --threads {threads} \
            --input {input.mags} \
            --extension .fa \
            --output-directory {params.out_dir} \
            --database_path {input.db}/uniref100.KO.1.dmnd \
            --remove_intermediates \
        2>> {log} 1>&2

        mv {params.out_dir}/quality_report.tsv {output} 2>> {log} 1>&2
        rm --recursive --verbose --force {params.out_dir} 2>> {log} 1>&2
        """


rule dereplicate__checkm2:
    """Run CheckM2"""
    input:
        rules._dereplicate__checkm2__predict.output,
