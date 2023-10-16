rule assemble_eval_cram_to_bam_one:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
        crai=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram.crai",
        reference=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        bam=temp(ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam"),
    log:
        ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam.log",
    conda:
        "assemble.yml"
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


rule assemble_eval_cram_to_bam_all:
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule assemble_eval_coverm_contig_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        tsv=ASSEMBLE_COVERM
        / "contig/{method}/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "assemble.yml"
    log:
        ASSEMBLE_COVERM / "contig/{method}/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["assemble"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["assemble"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.method} \
            --proper-pairs-only \
        > {output.tsv} 2> {log}
        """


rule assemble_eval_coverm_aggregate_contig:
    input:
        get_tsvs_for_assembly_coverm_genome,
    output:
        tsv=ASSEMBLE_COVERM / "contig.{method}.tsv",
    log:
        ASSEMBLE_COVERM / "contig.{method}.log",
    conda:
        "assemble.yml"
    params:
        input_dir=ASSEMBLE_COVERM / "contig/{method}",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule assemble_eval_coverm_contig:
    input:
        [
            ASSEMBLE_COVERM / f"contig.{method}.tsv"
            for method in params["assemble"]["coverm"]["contig"]["methods"]
        ],


rule assemble_eval_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        tsv=ASSEMBLE_COVERM
        / "genome/{method}/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "assemble.yml"
    log:
        ASSEMBLE_COVERM / "genome/{method}/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["assemble"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["assemble"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.method} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} 2> {log}
        """


rule assemble_eval_coverm_aggregate_genome:
    input:
        get_tsvs_for_assembly_coverm_genome,
    output:
        tsv=ASSEMBLE_COVERM / "genome.{method}.tsv",
    log:
        ASSEMBLE_COVERM / "genome.{method}.log",
    conda:
        "assemble.yml"
    params:
        input_dir=ASSEMBLE_COVERM / "genome",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule assemble_eval_coverm_genome:
    input:
        [
            ASSEMBLE_COVERM / f"genome.{method}.tsv"
            for method in params["assemble"]["coverm"]["genome"]["methods"]
        ],


rule assemble_eval_quast_one:
    """Run quast over one assembly group"""
    input:
        ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        directory(ASSEMBLE_QUAST / "{assembly_id}"),
    log:
        ASSEMBLE_QUAST / "{assembly_id}.log",
    conda:
        "assemble.yml"
    threads: 4
    params:
        extra=params["assemble"]["quast"]["extra"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """


rule assemble_eval_quast:
    """Run quast over all assembly groups"""
    input:
        [ASSEMBLE_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule assemble_eval_samtools:
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.{extension}"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
        ],


rule assemble_eval:
    input:
        rules.assemble_eval_quast.input,
        rules.assemble_eval_coverm_contig.input,
        rules.assemble_eval_coverm_genome.input,
        rules.assemble_eval_samtools.input,
