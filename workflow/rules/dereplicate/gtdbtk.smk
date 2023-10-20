rule dereplicate_eval_gtdbtk:
    input:
        bin_folder=DREP / "dereplicated_genomes",
        database=features["gtdbtk_database"],
    output:
        outdir=directory(DREP_GTDBTK),
    log:
        DREP_GTDBTK / "gtdbtk.log",
    conda:
        "gtdbtk.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["dereplicate"]["gtdbtk"]["memory_gb"]),
        runtime=24 * 60,
    shell:
        """
        export GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {input.bin_folder} \
            --extension fa \
            --out_dir {output.outdir} \
            --cpus {threads} \
            --skip_ani_screen \
            --full_tree \
        2> {log} 1>&2
        """


# coverm ----
rule dereplicate_eval_cram_to_bam_one:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa",
    output:
        bam=temp(DREP_BOWTIE2 / "{sample_id}.{library_id}.bam"),
    log:
        DREP_BOWTIE2 / "{sample_id}.{library_id}.bam.log",
    conda:
        "dereplicate.yml"
    resources:
        runtime=1 * 60,
        mem_mb=4 * 1024,
    shell:
        """
        samtools view \
            -F 4 \
            --reference {input.reference} \
            -o {output.bam} \
            -1 \
            {input.cram} \
        2> {log}
        """
