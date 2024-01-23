rule _checkv_run:
    input:
        fna=GENOMAD / "{assembly_id}_summary" / "{assembly_id}_virus.fna",
        database=features["databases"]["checkv"],
    output:
        fna=CHECKV / "{assembly_id}" / "combined.fna",
        summary=CHECKV / "{assembly_id}" / "quality_summary.tsv",
    log:
        CHECKV / "{assembly_id}/run.log"
    conda:
        "__environment__.yml"
    params:
        workdir=lambda w: CHECKV / f"{w.assembly_id}",
    threads:
        8
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


rule _checkv_filter:
    input:
        fna=CHECKV / "{assembly_id}" / "combined.fna",
        summary=CHECKV / "{assembly_id}" / "quality_summary.tsv",
    output:
        fna=CHECKV / "{assembly_id}" / "filtered.fna",
    log:
        CHECKV / "{assembly_id}" / "filter.log"
    params:
        min_quality=params["viral"]["checkv"]["min_quality"],
    conda:
        "__environment__.yml"
    shell:
        """
        medium='--regexp="Medium-quality"'
        high='--regexp="High-quality"'
        complete='--regexp="Complete"'
        
        if [[ "{params.min_quality}" == "Medium-quality" ]] ; then
            eval grep $medium $high $complete {input.summary}
        elif [[ "{params.min_quality}" == "High-quality" ]] ; then
            eval grep $high $complete {input.summary}
        else
            eval grep $complete {input.summary}
        fi \
        | cut \
            --fields 1 \
        | seqtk subseq \
            {input.fna} \
            /dev/stdin \
        > {output.fna}
        """


rule viral__checkv:
    input:
        [CHECKV / f"{assembly_id}" / "filtered.fna" for assembly_id in ASSEMBLIES]


