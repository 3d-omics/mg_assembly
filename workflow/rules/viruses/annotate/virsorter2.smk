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
    # shadow:
    #     "minimal"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    params:
        workdir=lambda w: VIR_VIRSORTER2 / w.assembly_id,
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


use rule csvtk__concat as viruses__annotate__virsorter2__concatenate_viral_boundary with:
    input:
        [
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-boundary.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_VIRSORTER2 / "final-viral-boundary.tsv.gz",
    log:
        VIR_VIRSORTER2 / "final-viral-boundary.log",


use rule csvtk__concat as viruses__annotate__virsorter2__concatenate_viral_score with:
    input:
        [
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-score.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_VIRSORTER2 / "final-viral-score.tsv.gz",
    log:
        VIR_VIRSORTER2 / "final-viral-score.log",


use rule csvtk__concat as viruses__annotate__virsorter2__concatenate_viral_contigs_for_dramv with:
    input:
        [
            VIR_VIRSORTER2 / f"{assembly_id}" / "viral-affi-contigs-for-dramv.tab"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_VIRSORTER2 / "viral-affi-contigs-for-dramv.tab.gz",
    log:
        VIR_VIRSORTER2 / "viral-affi-contigs-for-dramv.log",


use rule concatenate__gzip_text_files as viruses__annotate__virsorter2__concatenate_viral_combined with:
    input:
        [
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-combined.fa"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_VIRSORTER2 / "final-viral-combined.fa.gz",
    log:
        VIR_VIRSORTER2 / "final-viral-combined.log",


use rule concatenate__gzip_text_files as viruses__annotate__virsorter2__concatenate_viral_combined_for_dramv with:
    input:
        [
            VIR_VIRSORTER2 / f"{assembly_id}" / "final-viral-combined-for-dramv.fa"
            for assembly_id in ASSEMBLIES
        ]
        + NULL,
    output:
        VIR_VIRSORTER2 / "final-viral-combined-for-dramv.fa.gz",
    log:
        VIR_VIRSORTER2 / "final-viral-combined-for-dramv.log",


rule viruses__annotate__virsorter2__all:
    input:
        VIR_VIRSORTER2 / "final-viral-boundary.tsv.gz",
        VIR_VIRSORTER2 / "final-viral-score.tsv.gz",
        VIR_VIRSORTER2 / "viral-affi-contigs-for-dramv.tab.gz",
        VIR_VIRSORTER2 / "final-viral-combined.fa.gz",
        VIR_VIRSORTER2 / "final-viral-combined-for-dramv.fa.gz",
