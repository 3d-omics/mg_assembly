rule viruses__annotate__virsorter2__download:
    output:
        directory(features["databases"]["virsorter2"]),
    log:
        f"{features["databases"]["virsorter2"]}.log",
    conda:
        ENVS / "virsorter2.yml"
    shell:
        """
        virsorter setup \
            --db-dir {output} \
            --skip-deps-install \
            --jobs {threads} \
        2> {log} 1>&2
        """


rule viruses__annotate__virsorter2__run:
    input:
        fasta=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
        database=features["databases"]["virsorter2"],
    output:
        viruses_boundary=VIR_VIRSORTER2 / "{assembly_id}" / "final-viral-boundary.tsv",
        combined=VIR_VIRSORTER2 / "{assembly_id}" / "final-viral-combined.fa",
        score=VIR_VIRSORTER2 / "{assembly_id}" / "final-viral-score.tsv",
        dram_fa=VIR_VIRSORTER2 / "{assembly_id}" / "final-viral-combined-for-dramv.fa",
        dram_tsv=VIR_VIRSORTER2 / "{assembly_id}" / "viral-affi-contigs-for-dramv.tab",
    log:
        VIR_VIRSORTER2 / "{assembly_id}.log",
    conda:
        ENVS / "virsorter2.yml"
    params:
        workdir=lambda w: VIR_VIRSORTER2 / w.assembly_id,
    # shadow:
    #     "minimal"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    shell:
        """
        virsorter run \
            --working-dir {params.workdir} \
            --jobs {threads} \
            --prep-for-dramv \
            --tmpdir {params.workdir}/tmp \
            --rm-tmpdir \
            --verbose \
            --use-conda-off \
            --seqfile {input.fasta} \
            --db-dir {input.database} \
        2> {log} 1>&2

        mv \
            {params.workdir}/for-dramv/viral-affi-contigs-for-dramv.tab \
            {params.workdir}/for-dramv/final-viral-combined-for-dramv.fa \
            {params.workdir}/ \
        2>> {log} 1>&2

        rm -rfv {params.workdir}/for-dramv
        """


rule viruses__annotate__virsorter2__run__all:
    input:
        [
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-boundary.tsv"
            for assembly_id in ASSEMBLIES
        ],


rule viruses__annotate__virsorter2__aggregate_tsvs:
    input:
        boundary=[
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-boundary.tsv"
            for assembly_id in ASSEMBLIES
        ],
        score=[
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-score.tsv"
            for assembly_id in ASSEMBLIES
        ],
        contigs=[
            VIR_VIRSORTER2 / f"{assembly_id}" / "viral-affi-contigs-for-dramv.tab"
            for assembly_id in ASSEMBLIES
        ],
    output:
        boundary=VIR_VIRSORTER2 / "final-viral-boundary.tsv.gz",
        score=VIR_VIRSORTER2 / "final-viral-score.tsv.gz",
        contigs=VIR_VIRSORTER2 / "viral-affi-contigs-for-dramv.tab.gz",
    log:
        VIR_VIRSORTER2 / "aggregate_tsvs.log",
    conda:
        ENVS / "virsorter2.yml"
    threads: 24
    shell:
        """
        (
            csvtk concat \
                --tabs \
                {input.boundary} \
            | bgzip --compress-level 0 --threads {threads} \
            > {output.boundary}

            csvtk concat \
                --tabs \
                {input.score} \
            | bgzip --compress-level 0 --threads {threads} \
            > {output.score}

            csvtk concat \
                --tabs \
                {input.contigs} \
            | bgzip --compress-level 0 --threads {threads} \
            > {output.contigs}
        ) 2> {log}
        """


rule viruses__annotate__virsorter2__concatenate_fastas:
    input:
        combined=[
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-combined.fa"
            for assembly_id in ASSEMBLIES
        ],
        dramv_fa=[
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-combined-for-dramv.fa"
            for assembly_id in ASSEMBLIES
        ],
    output:
        combined=VIR_VIRSORTER2 / "final-viral-combined.fa.gz",
        dramv_fa=VIR_VIRSORTER2 / "final-viral-combined-for-dramv.fa.gz",
    log:
        VIR_VIRSORTER2 / "concatenate_fastas.log",
    conda:
        ENVS / "virsorter2.yml"
    shell:
        """
        (
            cat {input.combined} | bgzip --threads {threads} > {output.combined} 2> {log}
            cat {input.dramv_fa} | bgzip --threads {threads} > {output.dramv_fa} 2> {log}
        ) 2> {log}
        """


rule viruses__annotate__virsorter2__all:
    input:
        rules.viruses__annotate__virsorter2__aggregate_tsvs.output,
        rules.viruses__annotate__virsorter2__concatenate_fastas.output,
