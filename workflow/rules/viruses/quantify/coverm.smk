rule viruses__quantify__coverm__genome:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=VBOWTIE2 / "{sample_id}.{library_id}.bam",
        bai=VBOWTIE2 / "{sample_id}.{library_id}.bam.bai",
    output:
        tsv=VCOVERM / "genome" / "{method}" / "{sample_id}.{library_id}.tsv.gz",
    conda:
        "__environment__.yml"
    log:
        VCOVERM / "genome" / "{method}" / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["quantify"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["quantify"]["coverm"]["genome"]["separator"],
    shell:
        """
        ( coverm genome \
            --bam-files {input.bam} \
            --methods {params.method} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        | bgzip \
        > {output.tsv} \
        ) 2> {log}
        """


rule viruses__quantify__coverm__genome_aggregate:
    """Run coverm genome and a single method"""
    input:
        get_tsvs_for_dereplicate_vcoverm_genome,
    output:
        tsv=VCOVERM / "genome.{method}.tsv.gz",
    log:
        VCOVERM / "genome.{method}.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=lambda w: VCOVERM / "genome" / w.method,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule viruses__quantify__coverm__genome__all:
    """Run coverm genome and all methods"""
    input:
        [
            VCOVERM / f"genome.{method}.tsv.gz"
            for method in params["quantify"]["coverm"]["genome"]["methods"]
        ],


# coverm contig ----
rule viruses__quantify__coverm__contig:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=VBOWTIE2 / "{sample_id}.{library_id}.bam",
        bai=VBOWTIE2 / "{sample_id}.{library_id}.bam.bai",
    output:
        tsv=VCOVERM / "contig" / "{method}" / "{sample_id}.{library_id}.tsv.gz",
    conda:
        "__environment__.yml"
    log:
        VCOVERM / "contig" / "{method}" / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
    shell:
        """
        ( coverm contig \
            --bam-files {input.bam} \
            --methods {params.method} \
            --proper-pairs-only \
        | bgzip \
        > {output.tsv} \
        ) 2> {log}
        """


rule viruses__quantify__coverm__contig_aggregate:
    """Run coverm contig and a single method"""
    input:
        get_tsvs_for_dereplicate_vcoverm_contig,
    output:
        tsv=VCOVERM / "contig.{method}.tsv.gz",
    log:
        VCOVERM / "contig.{method}.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=lambda w: VCOVERM / "contig" / w.method,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule viruses__quantify__coverm__contig__all:
    """Run coverm contig and all methods"""
    input:
        [
            VCOVERM / f"contig.{method}.tsv.gz"
            for method in params["quantify"]["coverm"]["contig"]["methods"]
        ],


rule viruses__quantify__coverm__all:
    input:
        rules.viruses__quantify__coverm__contig__all.input,
        rules.viruses__quantify__coverm__genome__all.input,
