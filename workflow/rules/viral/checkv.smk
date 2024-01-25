rule _checkv_run:
    input:
        fna=GENOMAD / "{assembly_id}_summary" / "{assembly_id}_virus.fna",
        database=features["databases"]["checkv"],
    output:
        fna=CHECKV / "{assembly_id}" / "all.fna",
        summary=CHECKV / "{assembly_id}" / "quality_summary.tsv",
    log:
        CHECKV / "{assembly_id}/run.log",
    conda:
        "__environment__.yml"
    params:
        workdir=lambda w: CHECKV / f"{w.assembly_id}",
    threads: 8
    shell:
        """
        checkv end_to_end \
            {input.fna} \
            {params.workdir} \
            --remove_tmp \
            -t {threads} \
            --restart \
            -d {input.database} \
        2> {log} 1>&2

        cat \
            {params.workdir}/proviruses.fna \
            {params.workdir}/viruses.fna \
        > {output.fna} \
        2>> {log}
        """


rule viral__checkv:
    input:
        [CHECKV / f"{assembly_id}" / "all.fna" for assembly_id in ASSEMBLIES],
