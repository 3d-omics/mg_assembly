rule viruses__annotate__checkv__download:
    output:
        features["databases"]["checkv"],
    log:
        f"{features["databases"]["checkv"]}.log",
    shell:
        """
        checkv download_database \
            $(dirname {output}) \
        2> {log} 1>&2

        mv \
            --verbose \
            $(dirname {output}/checkv-db-v1.5) \
            {output} \
        2>> {log} 1>&2
        """


rule viruses__annotate__checkv__end_to_end:
    input:
        fasta=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
        database=features["databases"]["checkv"],
    output:
        complete_genomes=temp(VIR_CHECKV / "{assembly_id}" / "complete_genomes.tsv"),
        completeness=temp(VIR_CHECKV / "{assembly_id}" / "completeness.tsv"),
        contamination=temp(VIR_CHECKV / "{assembly_id}" / "contamination.tsv"),
        proviruses=temp(VIR_CHECKV / "{assembly_id}" / "proviruses.fna"),
        summary=temp(VIR_CHECKV / "{assembly_id}" / "quality_summary.tsv"),
        viruses=temp(VIR_CHECKV / "{assembly_id}" / "viruses.fna"),
    log:
        VIR_CHECKV / "{assembly_id}" / "checkv.log",
    conda:
        ENVS / "checkv.yml"
    params:
        workdir=lambda w: VIR_CHECKV / f"{w.assembly_id}",
    threads: 24
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    shell:
        """
        checkv end_to_end \
            -d {input.database} \
            -t {threads} \
            --restart \
            {input.fasta} \
            {params.workdir} \
        2> {log} 1>&2
        """


rule viruses__annotate__checkv__end_to_end__all:
    input:
        [
            VIR_CHECKV / f"{assembly_id}" / f"{file}"
            for assembly_id in ASSEMBLIES
            for file in [
                "complete_genomes.tsv",
                "completeness.tsv",
                "contamination.tsv",
                "quality_summary.tsv",
                "proviruses.fna",
                "viruses.fna",
            ]
        ],


rule viruses__annotate__checkv__aggregate_tsvs:
    input:
        complete_genomes=[
            VIR_CHECKV / f"{assembly_id}" / "complete_genomes.tsv"
            for assembly_id in ASSEMBLIES
        ],
        completeness=[
            VIR_CHECKV / f"{assembly_id}" / "completeness.tsv"
            for assembly_id in ASSEMBLIES
        ],
        contamination=[
            VIR_CHECKV / f"{assembly_id}" / "contamination.tsv"
            for assembly_id in ASSEMBLIES
        ],
        summary=[
            VIR_CHECKV / f"{assembly_id}" / "quality_summary.tsv"
            for assembly_id in ASSEMBLIES
        ],
    output:
        complete_genomes=VIR_CHECKV / "checkv.complete_genomes.tsv.gz",
        completeness=VIR_CHECKV / "checkv.completeness.tsv.gz",
        contamination=VIR_CHECKV / "checkv.contamination.tsv.gz",
        summary=VIR_CHECKV / "checkv.quality_summary.tsv.gz",
    log:
        VIR_CHECKV / "checkv.aggregate_tsvs.log",
    conda:
        ENVS / "checkv.yml"
    shell:
        """
        ( csvtk concat --tabs {input.complete_genomes} \
        | bgzip --compress-level 9 --threads {threads} \
        > {output.complete_genomes} ) 2> {log}

        ( csvtk concat --tabs {input.completeness} \
        | bgzip --compress-level 9 --threads {threads} \
        > {output.completeness} ) 2>> {log}

        ( csvtk concat --tabs {input.contamination} \
        | bgzip --compress-level 9 --threads {threads} \
        > {output.contamination} ) 2>> {log}

        ( csvtk concat --tabs {input.summary} \
        | bgzip --compress-level 9 --threads {threads} \
        > {output.summary} ) 2>> {log}
        """


rule viruses__annotate__checkv__concatenate_fastas:
    input:
        proviruses=[
            VIR_CHECKV / f"{assembly_id}" / "proviruses.fna"
            for assembly_id in ASSEMBLIES
        ],
        viruses=[
            VIR_CHECKV / f"{assembly_id}" / "viruses.fna" for assembly_id in ASSEMBLIES
        ],
    output:
        proviruses=VIR_CHECKV / "checkv.proviruses.fna.gz",
        viruses=VIR_CHECKV / "checkv.viruses.fna.gz",
    log:
        VIR_CHECKV / "checkv.concatenate_fastas.log",
    conda:
        ENVS / "checkv.yml"
    shell:
        """
        ( cat {input.proviruses} \
        | bgzip --compress-level 9 --threads {threads} \
        > {output.proviruses} \
        ) 2> {log}

        ( cat {input.viruses} \
        | bgzip --compress-level 9 --threads {threads} \
        > {output.viruses} \
        ) 2>> {log}
        """


rule viruses__annotate__checkv__all:
    input:
        rules.viruses__annotate__checkv__aggregate_tsvs.output,
        rules.viruses__annotate__checkv__concatenate_fastas.output,
