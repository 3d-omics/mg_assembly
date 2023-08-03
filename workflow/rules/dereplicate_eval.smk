
rule dereplicate_eval_gtdbtk:
    input:
        bin_folder=DREP / "dereplicated_genomes",
        database=features["gtdbtk_database"],
    output:
        outdir=DREP_GTDBTK,
        bac120=DREP_GTDBTK / "classify/gtdbtk.bac120.summary.tsv",
        combined_summary=DREP_GTDBTK / "_gtdbtk_combined_summary.tsv",
        tree=DREP_GTDBTK / "_gtdbtk.bac120.classify.tree",
        taxonomy="_taxon_table.tsv",
    log:
        DREP_GTDBTK / "gtdbtk.log",
    conda:
        "../envs/dereplicate.yml"
    threads: 24
    resources:
        mem_mb=150 * 1024,
    shell:
        """
        GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {params.bins} \
            --extension gz \
            --out_dir {output.outdir} \
            --cpus {threads} \
            --skip_ani_screen \
            --full_tree \
        2> {log} 1>&2
        """


rule dereplicate_eval_cram_to_bam_one:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=DREP_BOWTIE2 / "dereplicated_genomes.{sample_id}.{library_id}.cram",
        crai=DREP_BOWTIE2 / "dereplicated_genomes.{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa",
    output:
        bam=temp(DREP_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam"),
    log:
        DREP_BOWTIE2 / "{assembly_id},{sample_id}.{library_id}.bam.log",
    conda:
        "../envs/dereplicate.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=4 * 1024,
    shell:
        """
        samtools view \
            -F 4 \
            --threads {threads} \
            --reference {input.reference} \
            --output {output.bam} \
            --fast \
            {input.cram} \
        2> {log}
        """


rule dereplicate_eval_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=DREP_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
    output:
        tsv=DREP_COVERM / "genome/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/dereplicate.yml"
    log:
        DREP_COVERM / "genome/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["dereplicate"]["coverm"]["genome"]["methods"],
        min_covered_fraction=params["dereplicate"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["dereplicate"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} 2> {log}
        """


rule dereplicate_eval_coverm_genome:
    input:
        tsvs=[
            DREP_COVERM / f"genome/{assembly_id}.{sample_id}.{library_id}.tsv"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
    output:
        tsv=DREP_COVERM / "genome.tsv",
    log:
        DREP_COVERM / "genome.log",
    conda:
        "../envs/assembly.yml"
    params:
        input_dir=DREP_COVERM / "genome/",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """
