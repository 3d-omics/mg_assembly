rule prokaryotes__quantify__coverm__genome__:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        cram=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa.gz",
        fai=DREP / "dereplicated_genomes.fa.gz.fai",
    output:
        tsv=COVERM / "genome" / "{method}" / "{sample_id}.{library_id}.tsv.gz",
    conda:
        "__environment__.yml"
    log:
        COVERM / "genome" / "{method}" / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["quantify"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["quantify"]["coverm"]["genome"]["separator"],
    shell:
        """
        ( samtools view \
            --exclude-flags 4 \
            --reference {input.reference} \
            --fast \
            {input.cram} \
        | coverm genome \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        | gzip \
        > {output.tsv} \
        ) 2> {log}
        """


rule prokaryotes__quantify__coverm__genome__aggregate__:
    """Run coverm genome and a single method"""
    input:
        get_tsvs_for_dereplicate_coverm_genome,
    output:
        tsv=COVERM / "genome.{method}.tsv.gz",
    log:
        COVERM / "genome.{method}.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=lambda w: COVERM / "genome" / w.method,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule prokaryotes__quantify__coverm__genome:
    """Run coverm genome and all methods"""
    input:
        [
            COVERM / f"genome.{method}.tsv.gz"
            for method in params["quantify"]["coverm"]["genome"]["methods"]
        ],


# coverm contig ----
rule prokaryotes__quantify__coverm__contig__:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        cram=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa.gz",
        fai=DREP / "dereplicated_genomes.fa.gz.fai",
    output:
        tsv=COVERM / "contig" / "{method}" / "{sample_id}.{library_id}.tsv.gz",
    conda:
        "__environment__.yml"
    log:
        COVERM / "contig" / "{method}" / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
    shell:
        """
        ( samtools view \
            --exclude-flags 4 \
            --reference {input.reference} \
            --fast \
            {input.cram} \
        | coverm contig \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --proper-pairs-only \
        | gzip \
        > {output.tsv} \
        ) 2> {log}
        """


rule prokaryotes__quantify__coverm__contig__aggregate__:
    """Run coverm contig and a single method"""
    input:
        get_tsvs_for_dereplicate_coverm_contig,
    output:
        tsv=COVERM / "contig.{method}.tsv.gz",
    log:
        COVERM / "contig.{method}.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=lambda w: COVERM / "contig" / w.method,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule prokaryotes__quantify__coverm__contig:
    """Run coverm contig and all methods"""
    input:
        [
            COVERM / f"contig.{method}.tsv.gz"
            for method in params["quantify"]["coverm"]["contig"]["methods"]
        ],


rule prokaryotes__quantify__coverm:
    input:
        rules.prokaryotes__quantify__coverm__genome.input,
        rules.prokaryotes__quantify__coverm__contig.input,
