# gtdbtk ----
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
        mem_mb=320 * 1024,
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


## coverm genome ----
rule dereplicate_eval_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=DREP_BOWTIE2 / "{sample_id}.{library_id}.bam",
    output:
        tsv=DREP_COVERM / "genome/{sample_id}.{library_id}.tsv",
    conda:
        "dereplicate.yml"
    log:
        DREP_COVERM / "genome/{sample_id}.{library_id}.log",
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
            DREP_COVERM / f"genome/{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        tsv=DREP_COVERM / "genome.tsv",
    log:
        DREP_COVERM / "genome.log",
    conda:
        "dereplicate.yml"
    params:
        input_dir=DREP_COVERM / "genome/",
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


# coverm contig ----
rule dereplicate_eval_coverm_contig_one:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=DREP_BOWTIE2 / "{sample_id}.{library_id}.bam",
    output:
        tsv=DREP_COVERM / "contig/{sample_id}.{library_id}.tsv",
    conda:
        "dereplicate.yml"
    log:
        DREP_COVERM / "contig/{sample_id}.{library_id}.log",
    params:
        methods=params["dereplicate"]["coverm"]["contig"]["methods"],
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output.tsv} 2> {log}
        """


rule dereplicate_eval_coverm_contig:
    input:
        tsvs=[
            DREP_COVERM / f"contig/{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        tsv=DREP_COVERM / "contig.tsv",
    log:
        DREP_COVERM / "contig.log",
    conda:
        "dereplicate.yml"
    params:
        input_dir=DREP_COVERM / "contig/",
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


# dram ----
rule dereplicate_eval_dram_annotate:
    input:
        drep_folder=DREP / "dereplicated_genomes",
        mock_db="results/dram_db_setup.done",
    output:
        outdir=DREP_DRAM / "annotate/dereplicated_genomes",
        annotations=DREP_DRAM / "annotate/dereplicated_genomes/annotations.tsv",
        trnas=touch(DREP_DRAM / "annotate/dereplicated_genomes/trnas.tsv"),
        rrnas=touch(DREP_DRAM / "annotate/dereplicated_genomes/rrnas.tsv"),
    log:
        DREP_DRAM / "annotate/dereplicated_genomes.log",
    conda:
        "dram.yml"
    params:
        min_contig_size=1500,
    threads: 24
    shell:
        """
        rm -rf {output.outdir}

        DRAM.py annotate \
            --input_fasta {input.drep_folder} \
            --output_dir {output.outdir} \
            --threads {threads} \
            --rrna_path {output.rrnas} \
            --trna_path {output.trnas} \
            --min_contig_size {params.min_contig_size} \
        2> {log} 1>&2
        """


rule dereplicate_eval_dram_distill:
    input:
        indir=DREP_DRAM / "annotate/dereplicated_genomes",
        annotations=DREP_DRAM / "annotate/dereplicated_genomes/annotations.tsv",
        trnas=DREP_DRAM / "annotate/dereplicated_genomes/trnas.tsv",
        rrnas=DREP_DRAM / "annotate/dereplicated_genomes/rrnas.tsv",
        mock_db="results/dram_db_setup.done",
    output:
        outdir=DREP_DRAM / "distill/dereplicated_genomes",
    log:
        DREP_DRAM / "distill/dereplicated_genomes.log",
    conda:
        "dram.yml"
    shell:
        """
        DRAM.py distill \
            --input_dir {input.indir} \
            --output_dir {output.outdir} \
            --annotations {input.annotations} \
            --rrna_path {input.rrnas} \
            --trna_path {input.trnas} \
        2> {log} 1>&2
        """


# quast ----
rule dereplicate_eval_quast:
    """Run quast over one the dereplicated mags"""
    input:
        DREP / "dereplicated_genomes.fa",
    output:
        directory(DREP_QUAST / "dereplicated_genomes"),
    log:
        DREP_QUAST / "dereplicated_genomes.log",
    conda:
        "dereplicate.yml"
    threads: 4
    params:
        extra=params["metabin"]["quast"]["extra"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """


# samtools ----
rule dereplicate_eval_samtools:
    input:
        [
            DREP_BOWTIE2 / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
        ],


# eval ----
rule dereplicate_eval:
    input:
        rules.dereplicate_eval_coverm_genome.output,
        rules.dereplicate_eval_coverm_contig.output,
        rules.dereplicate_eval_quast.output,


rule dereplicate_eval_with_dram:
    input:
        rules.dereplicate_eval.input,
        rules.dereplicate_eval_dram_distill.output,


rule dereplicate_eval_with_gtdbtk:
    input:
        rules.dereplicate_eval.input,
        rules.dereplicate_eval_gtdbtk.output,
