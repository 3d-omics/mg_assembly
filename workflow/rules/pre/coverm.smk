rule pre_cram_to_mapped_bam:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=get_cram_for_pre_eval_cram_to_mapped_bam,
        crai=get_crai_for_pre_eval_cram_to_mapped_bam,
        reference=REFERENCE / f"{LAST_HOST}.fa.gz",
        fai=REFERENCE / f"{LAST_HOST}.fa.gz.fai",
    output:
        bam=temp(PRE_COVERM / "bams" / "{sample_id}.{library_id}.bam"),
    log:
        PRE_COVERM / "bams" / "{sample_id}.{library_id}.log",
    conda:
        "_env.yml"
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


rule pre_coverm_genome_method_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=PRE_COVERM / "bams" / "{sample}.{library_id}.bam",
        bai=PRE_COVERM / "bams" / "{sample}.{library_id}.bam.bai",
    output:
        tsv=touch(PRE_COVERM / "genome" / "{method}" / "{sample}.{library_id}.tsv"),
    conda:
        "_env.yml"
    log:
        PRE_COVERM / "genome" / "{method}" / "{sample}.{library_id}.log",
    params:
        methods="{method}",
        min_covered_fraction=params["pre"]["coverm"]["genome"]["min_covered_fraction"],
        separator=params["pre"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output} 2> {log} || true \
        """


rule pre_coverm_genome_aggregate_method:
    """Aggregate all the nonpareil results into a single table"""
    input:
        [
            PRE_COVERM / "genome" / "{method}" / f"{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        PRE_COVERM / "genome.{method}.tsv",
    log:
        PRE_COVERM / "genome.{method}.log",
    conda:
        "_env.yml"
    params:
        input_dir=PRE_COVERM / "genome" / "{method}",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule pre_coverm:
    """Run both coverm genome over the preprocessed samples"""
    input:
        [
            PRE_COVERM / f"genome.{method}.tsv"
            for method in params["pre"]["coverm"]["genome"]["methods"]
        ],
