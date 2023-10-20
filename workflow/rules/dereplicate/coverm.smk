rule dereplicate_eval_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=DREP_BOWTIE2 / "{sample_id}.{library_id}.bam",
    output:
        tsv=DREP_COVERM / "genome/{method}/{sample_id}.{library_id}.tsv",
    conda:
        "dereplicate.yml"
    log:
        DREP_COVERM / "genome/{method}/{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["dereplicate"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["dereplicate"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.method} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} 2> {log}
        """


rule dereplicate_eval_coverm_genome_method:
    input:
        get_tsvs_for_dereplicate_coverm_genome,
    output:
        tsv=DREP_COVERM / "genome.{method}.tsv",
    log:
        DREP_COVERM / "genome.{method}.log",
    conda:
        "dereplicate.yml"
    params:
        input_dir=DREP_COVERM / "genome/{method}",
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule dereplicate_eval_coverm_genome:
    input:
        [
            DREP_COVERM / f"genome.{method}.tsv"
            for method in params["dereplicate"]["coverm"]["genome"]["methods"]
        ],


# coverm contig ----
rule dereplicate_eval_coverm_contig_method_one:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=DREP_BOWTIE2 / "{sample_id}.{library_id}.bam",
    output:
        tsv=DREP_COVERM / "contig/{method}/{sample_id}.{library_id}.tsv",
    conda:
        "dereplicate.yml"
    log:
        DREP_COVERM / "contig/{method}/{sample_id}.{library_id}.log",
    params:
        method="{method}",
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.method} \
            --proper-pairs-only \
        > {output.tsv} 2> {log}
        """


rule dereplicate_eval_coverm_contig_method:
    input:
        get_tsvs_for_dereplicate_coverm_contig,
    output:
        tsv=DREP_COVERM / "contig.{method}.tsv",
    log:
        DREP_COVERM / "contig.{method}.log",
    conda:
        "dereplicate.yml"
    params:
        input_dir=DREP_COVERM / "contig/{method}",
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule dereplicate_eval_coverm_contig:
    input:
        [
            DREP_COVERM / f"contig.{method}.tsv"
            for method in params["dereplicate"]["coverm"]["contig"]["methods"]
        ],