rule viruses__annotate__checkv__download:
    output:
        features["databases"]["checkv"],
    log:
        f"{features["databases"]["checkv"]}.log",
    conda:
        "base"
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
    threads: 24
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    params:
        workdir=lambda w: VIR_CHECKV / f"{w.assembly_id}",
    shell:
        """
        checkv end_to_end \
            -d {input.database} \
            -t {threads} \
            --restart \
            {input.fasta} \
            {params.workdir} \
        2> {log} 1>&2

        rm \
            --recursive \
            --force \
            --verbose \
            {params.workdir}/tmp/ \
        2>> {log} 1>&2
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


use rule csvtk__concat as viruses__annotate__checkv__concatenate_complete_genomes with:
    input:
        [
            VIR_CHECKV / f"{assembly_id}" / "complete_genomes.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_CHECKV / "checkv.complete_genomes.tsv.gz",
    log:
        VIR_CHECKV / "checkv.complete_genomes.log",


use rule csvtk__concat as viruses__annotate__checkv__concatenate_completeness with:
    input:
        [
            VIR_CHECKV / f"{assembly_id}" / "completeness.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_CHECKV / "checkv.completeness.tsv.gz",
    log:
        VIR_CHECKV / "checkv.completeness.log",


use rule csvtk__concat as viruses__annotate__checkv__concatenate_contamination with:
    input:
        [
            VIR_CHECKV / f"{assembly_id}" / "contamination.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_CHECKV / "checkv.contamination.tsv.gz",
    log:
        VIR_CHECKV / "checkv.contamination.tsv.gz",


use rule csvtk__concat as viruses__annotate__checkv__concatenate_summary with:
    input:
        [
            VIR_CHECKV / f"{assembly_id}" / "quality_summary.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_CHECKV / "checkv.quality_summary.tsv.gz",
    log:
        VIR_CHECKV / "checkv.quality_summary.tsv.gz",


use rule concatenate__gzip_text_files as viruses__annotate__checkv__concatenate_proviruses with:
    input:
        [
            VIR_CHECKV / f"{assembly_id}" / "proviruses.fna"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_CHECKV / "checkv.proviruses.fna.gz",
    log:
        VIR_CHECKV / "checkv.proviruses.log",


use rule concatenate__gzip_text_files as viruses__annotate__checkv__concatenate_viruses with:
    input:
        [VIR_CHECKV / f"{assembly_id}" / "viruses.fna" for assembly_id in ASSEMBLIES],
    output:
        VIR_CHECKV / "checkv.viruses.fna.gz",
    log:
        VIR_CHECKV / "checkv.viruses.log",


rule viruses__annotate__checkv__all:
    input:
        rules.viruses__annotate__checkv__end_to_end__all.input,
        [
            VIR_CHECKV / "checkv.complete_genomes.tsv.gz",
            VIR_CHECKV / "checkv.completeness.tsv.gz",
            VIR_CHECKV / "checkv.contamination.tsv.gz",
            VIR_CHECKV / "checkv.quality_summary.tsv.gz",
            VIR_CHECKV / "checkv.proviruses.fna.gz",
            VIR_CHECKV / "checkv.viruses.fna.gz",
        ],
