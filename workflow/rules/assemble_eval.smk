rule assemble_eval_coverm_contig_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=ASSEMBLY_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=ASSEMBLY_RENAME / "{assembly_id}.fa",
    output:
        tsv=ASSEMBLY_COVERM / "contig/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/assembly.yml"
    log:
        ASSEMBLY_COVERM / "contig/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["assemble"]["coverm"]["genome"]["methods"],
        min_covered_fraction=params["assemble"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["assemble"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output.tsv} 2> {log}
        """


rule assemble_eval_coverm_aggregate_contig:
    input:
        tsvs=[
            ASSEMBLY_COVERM / f"contig/{assembly_id}.{sample_id}.{library_id}.tsv"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
    output:
        tsv=ASSEMBLY_COVERM / "contig.tsv",
    log:
        ASSEMBLY_COVERM / "contig.log",
    conda:
        "../envs/assembly.yml"
    params:
        input_dir=ASSEMBLY_COVERM / "contig",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule assemble_eval_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=ASSEMBLY_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=ASSEMBLY_RENAME / "{assembly_id}.fa",
    output:
        tsv=ASSEMBLY_COVERM / "genome/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/assembly.yml"
    log:
        ASSEMBLY_COVERM / "genome/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["assemble"]["coverm"]["genome"]["methods"],
        min_covered_fraction=params["assemble"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["assemble"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} 2> {log}
        """


rule assemble_eval_coverm_aggregate_genome:
    input:
        tsvs=[
            ASSEMBLY_COVERM / f"genome/{assembly_id}.{sample_id}.{library_id}.tsv"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
    output:
        tsv=ASSEMBLY_COVERM / "genome.tsv",
    log:
        ASSEMBLY_COVERM / "genome.log",
    conda:
        "../envs/assembly.yml"
    params:
        input_dir=ASSEMBLY_COVERM / "genome",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule assemble_eval_quast_one:
    """Run quast over one assembly group"""
    input:
        ASSEMBLY_RENAME / "{assembly_id}.fa",
    output:
        directory(ASSEMBLY_QUAST / "{assembly_id}"),
    log:
        ASSEMBLY_QUAST / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
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


rule assemble_eval_quast_all:
    """Run quast over all assembly groups"""
    input:
        [ASSEMBLY_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule assemble_eval_samtools:
    input:
        [
            ASSEMBLY_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.{extension}"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
        ],


rule assemble_eval:
    input:
        rules.assemble_eval_quast_all.input,
        rules.assemble_eval_coverm_aggregate_genome.output,
        rules.assemble_eval_coverm_aggregate_contig.output,
        rules.assemble_eval_samtools.input,
