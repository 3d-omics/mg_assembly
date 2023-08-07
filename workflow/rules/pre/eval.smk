rule pre_eval_nonpareil_one:
    """Run nonpareil over one sample

    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    Note3: in case that nonpareil fails for low coverage samples, it creates
    empty files
    """
    input:
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
    output:
        forward_fq=temp(NONPAREIL / "{sample_id}.{library_id}_1.fq"),
        npa=touch(NONPAREIL / "{sample_id}.{library_id}.npa"),
        npc=touch(NONPAREIL / "{sample_id}.{library_id}.npc"),
        npl=touch(NONPAREIL / "{sample_id}.{library_id}.npl"),
        npo=touch(NONPAREIL / "{sample_id}.{library_id}.npo"),
    log:
        NONPAREIL / "{sample_id}.{library_id}.log",
    conda:
        "pre.yml"
    params:
        prefix=compose_prefix_for_nonpareil,
    resources:
        runtime=24 * 60,
    shell:
        """
        gzip -dc {input.forward_} > {output.forward_fq} 2> {log}

        nonpareil \
            -s {output.forward_fq} \
            -T kmer \
            -b {params.prefix} \
            -f fastq \
            -t {threads} \
        2>> {log} 1>&2
        """


rule pre_eval_nonpareil:
    """Aggregate all the nonpareil results into a single table"""
    input:
        [
            NONPAREIL / f"{sample_id}.{library_id}.npo"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        NONPAREIL / "nonpareil.tsv",
    log:
        NONPAREIL / "nonpareil.log",
    conda:
        "pre.yml"
    params:
        input_dir=NONPAREIL,
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_nonpareil.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule pre_eval_singlem_pipe_one:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        data=features["singlem_database"],
    output:
        archive_otu_table=SINGLEM / "{sample_id}.{library_id}.archive.json",
        otu_table=SINGLEM / "{sample_id}.{library_id}.otu_table.tsv",
        condense=SINGLEM / "{sample_id}.{library_id}.condense.tsv",
    log:
        SINGLEM / "{sample_id}.{library_id}.log",
    conda:
        "pre.yml"
    threads: 4
    resources:
        runtime=4 * 60,
    shell:
        """
        singlem pipe \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --otu-table {output.otu_table} \
            --archive-otu-table {output.archive_otu_table} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
            --threads {threads} \
            --assignment-threads {threads} \
        2> {log} 1>&2 || true
        """


rule pre_eval_singlem:
    """Aggregate all the singlem results into a single table"""
    input:
        archive_otu_tables=[
            SINGLEM / f"{sample_id}.{library_id}.archive.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        database=features["singlem_database"],
    output:
        condense=SINGLEM / "singlem.tsv",
    log:
        SINGLEM / "singlem.log",
    conda:
        "pre.yml"
    params:
        input_dir=SINGLEM,
    shell:
        """
        singlem condense \
            --input-archive-otu-tables {input.archive_otu_tables} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.database} \
        2> {log} 1>&2
        """


rule pre_eval_cram_to_mapped_bam:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=PRE_BOWTIE2 / "{sample}.{library_id}.cram",
        reference=features["host"]["fasta"],
    output:
        bam=temp(PRE_COVERM / "bams/{sample}.{library_id}.bam"),
    log:
        PRE_COVERM / "bams/{sample}.{library_id}.log",
    conda:
        "pre.yml"
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


rule pre_eval_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=PRE_COVERM / "bams/{sample}.{library_id}.bam",
    output:
        tsv=touch(PRE_COVERM / "genome/{sample}.{library_id}.tsv"),
    conda:
        "pre.yml"
    log:
        PRE_COVERM / "genome/{sample}.{library_id}.log",
    params:
        methods=params["pre"]["coverm"]["genome"]["methods"],
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


rule pre_eval_coverm:
    """Aggregate all the nonpareil results into a single table"""
    input:
        [
            PRE_COVERM / f"genome/{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        PRE_COVERM / "coverm.tsv",
    log:
        PRE_COVERM / "coverm.log",
    conda:
        "pre.yml"
    params:
        input_dir=PRE_COVERM / "genome",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule pre_eval_samtools:
    input:
        [
            PRE_BOWTIE2 / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
        ],


rule pre_eval_fastqc_fastp:
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule pre_eval_fastqc_nonhost:
    input:
        [
            NONHOST / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule pre_eval:
    input:
        rules.pre_eval_nonpareil.output,
        rules.pre_eval_singlem.output,
        rules.pre_eval_coverm.output,
        rules.pre_eval_samtools.input,
        rules.pre_eval_fastqc_fastp.input,
        rules.pre_eval_fastqc_nonhost.input,
