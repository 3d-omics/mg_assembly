rule _dereplicate__coverm__cram_to_bam:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa",
        fai=DREP / "dereplicated_genomes.fa.fai",
    output:
        bam=temp(DREP_COVERM / "bams" / "{sample_id}.{library_id}.bam"),
    log:
        DREP_COVERM / "bams" / "{sample_id}.{library_id}.bam.log",
    conda:
        "_env.yml"
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


rule _dereplicate__coverm__genome:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=DREP_COVERM / "bams" / "{sample_id}.{library_id}.bam",
        bai=DREP_COVERM / "bams" / "{sample_id}.{library_id}.bam.bai",
    output:
        tsv=DREP_COVERM / "genome" / "{method}" / "{sample_id}.{library_id}.tsv",
    conda:
        "_env.yml"
    log:
        DREP_COVERM / "genome" / "{method}" / "{sample_id}.{library_id}.log",
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


rule _dereplicate__coverm__genome_aggregate:
    """Run coverm genome and a single method"""
    input:
        get_tsvs_for_dereplicate_coverm_genome,
    output:
        tsv=DREP_COVERM / "genome.{method}.tsv",
    log:
        DREP_COVERM / "genome.{method}.log",
    conda:
        "_env.yml"
    params:
        input_dir=compose_input_dir_for_dereplicate_coverm_genome_method,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule dereplicate__coverm__genome:
    """Run coverm genome and all methods"""
    input:
        [
            DREP_COVERM / f"genome.{method}.tsv"
            for method in params["dereplicate"]["coverm"]["genome"]["methods"]
        ],


# coverm contig ----
rule _dereplicate__coverm__contig:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=DREP_COVERM / "bams" / "{sample_id}.{library_id}.bam",
        bai=DREP_COVERM / "bams" / "{sample_id}.{library_id}.bam.bai",
    output:
        tsv=DREP_COVERM / "contig" / "{method}" / "{sample_id}.{library_id}.tsv",
    conda:
        "_env.yml"
    log:
        DREP_COVERM / "contig" / "{method}" / "{sample_id}.{library_id}.log",
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


rule _dereplicate__coverm__contig_aggregate:
    """Run coverm contig and a single method"""
    input:
        get_tsvs_for_dereplicate_coverm_contig,
    output:
        tsv=DREP_COVERM / "contig.{method}.tsv",
    log:
        DREP_COVERM / "contig.{method}.log",
    conda:
        "_env.yml"
    params:
        input_dir=compose_input_dir_for_dereplicate_coverm_contig_method,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule dereplicate__coverm__contig:
    """Run coverm contig and all methods"""
    input:
        [
            DREP_COVERM / f"contig.{method}.tsv"
            for method in params["dereplicate"]["coverm"]["contig"]["methods"]
        ],


rule dereplicate__coverm:
    input:
        rules.dereplicate__coverm__contig.input,
        rules.dereplicate__coverm__genome.input,
