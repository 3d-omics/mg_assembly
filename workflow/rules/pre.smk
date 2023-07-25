rule pre_fastp_trim_one:
    """Run fastp on one library"""
    input:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
    output:
        forward_=temp(FASTP / "{sample_id}.{library_id}_1.fq.gz"),
        reverse_=temp(FASTP / "{sample_id}.{library_id}_2.fq.gz"),
        unpaired1=temp(FASTP / "{sample_id}.{library_id}_u1.fq.gz"),
        unpaired2=temp(FASTP / "{sample_id}.{library_id}_u2.fq.gz"),
        html=FASTP / "{sample_id}.{library_id}_fastp.html",
        json=FASTP / "{sample_id}.{library_id}_fastp.json",
    log:
        FASTP / "{sample_id}.{library_id}.log",
    conda:
        "../envs/pre.yml"
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
        extra=params["fastp"]["extra"],
        length_required=params["fastp"]["length_required"],
    threads: 16
    resources:
        mem_mb=4 * 1024,
        runtime=240,
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 >(gzip --fast > {output.forward_}) \
            --out2 >(gzip --fast > {output.reverse_}) \
            --unpaired1 >(gzip --fast > {output.unpaired1}) \
            --unpaired2 >(gzip --fast > {output.unpaired2}) \
            --html {output.html} \
            --json {output.json} \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2
        """


rule pre_fastp_trim_all:
    """Run fastp over all libraries"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
        ],


rule pre_fastp_fastqc_one:
    """Run fastqc on one library from fastp output"""
    input:
        fq=FASTP / "{sample_id}.{library_id}_{end}.fq.gz",
    output:
        html=FASTP / "{sample_id}.{library_id}_{end}_fastqc.html",
        zip_=FASTP / "{sample_id}.{library_id}_{end}_fastqc.zip",
    log:
        FASTP / "{sample_id}.{library_id}_{end}_fastqc.log",
    conda:
        "../envs/pre.yml"
    shell:
        """
        fastqc \
            --outdir {FASTP} \
            --threads 1 \
            {input.fq} \
        2> {log} 1>&2
        """


rule pre_fastp_fastqc_all:
    """Run fastqc over all libraries after fastp"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
            for extension in "html zip".split(" ")
        ],


rule pre_bowtie2_index_host:
    """Build bowtie2 index for the human reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=features["host"]["fasta"],
    output:
        mock=touch(BOWTIE2 / "host"),
    log:
        BOWTIE2 / "host_index.log",
    conda:
        "../envs/pre.yml"
    params:
        extra=params["bowtie2"]["extra"],
    threads: 8
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.reference} \
            {output.mock} \
        2> {log} 1>&2
        """


rule pre_bowtie2_map_host_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        mock=BOWTIE2 / "host",
        reference=features["host"]["fasta"],
    output:
        cram=BOWTIE2 / "{sample_id}.{library_id}.cram",
    log:
        BOWTIE2 / "{sample_id}.{library_id}.log",
    conda:
        "../envs/pre.yml"
    params:
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        (bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule pre_bowtie2_map_host_all:
    """Map all libraries to reference genome using bowtie2"""
    input:
        [
            BOWTIE2 / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIB
        ],


# rule pre_bowtie2_reports


rule pre_bowtie2_extract_nonhost_one:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=BOWTIE2 / "{sample_id}.{library_id}.cram",
        reference=features["host"]["fasta"],
    output:
        forward_=temp(NONHOST / "{sample_id}.{library_id}_1.fq.gz"),
        reverse_=temp(NONHOST / "{sample_id}.{library_id}_2.fq.gz"),
    log:
        NONHOST / "{sample_id}.{library_id}.log",
    conda:
        "../envs/pre.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=32 * 1024,
    shell:
        """
        (samtools view \
            --reference {input.reference} \
            --threads {threads} \
            -u \
            -o /dev/stdout \
            -f 12 \
            {input.cram} \
        | samtools sort \
            -n \
            -u \
            --threads {threads} \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null \
            -c 1 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule pre_bowtie2_extract_nonhost_all:
    """Extract nonhost reads from all libraries"""
    input:
        [
            NONHOST / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIB
            for end in "1 2".split(" ")
        ],


rule pre_nonpareil_one:
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
        "../envs/pre.yml"
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


rule pre_nonpareil_all:
    """Run pre_nonpareil_one for all the samples"""
    input:
        [
            NONPAREIL / f"{sample_id}.{library_id}.{extension}"
            for extension in ["npa", "npc", "npl", "npo"]
            for sample_id, library_id in SAMPLE_LIB
        ],


rule pre_nonpareil:
    """Aggregate all the nonpareil results into a single table"""
    input:
        rules.pre_nonpareil_all.input,
    output:
        NONPAREIL / "nonpareil.tsv",
    log:
        NONPAREIL / "nonpareil.log",
    conda:
        "../envs/pre.yml"
    params:
        input_dir=NONPAREIL,
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_nonpareil.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule pre_singlem_data:
    """Download the singlem data

    For reasons unknown, you have to specify the filename, that may change in
    the future.
    """
    output:
        directory(SINGLEM / "data/S3.2.0.GTDB_r214.metapackage_20230428.smpkg.zb"),
    log:
        SINGLEM / "data.log",
    conda:
        "../envs/pre.yml"
    params:
        output_prefix=SINGLEM / "data",
    shell:
        """
        singlem data --output-directory {params.output_prefix} 2> {log} 1>&2
        """


rule pre_singlem_pipe_one:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        data=rules.pre_singlem_data.output,
    output:
        archive_otu_table=SINGLEM / "{sample_id}.{library_id}.archive.json",
        otu_table=SINGLEM / "{sample_id}.{library_id}.otu_table.tsv",
        condense=SINGLEM / "{sample_id}.{library_id}.condense.tsv",
    log:
        SINGLEM / "{sample_id}.{library_id}.log",
    conda:
        "../envs/pre.yml"
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


rule pre_singlem_pipe_all:
    """Run pre_singlem_one for all the samples"""
    input:
        [
            SINGLEM / f"{sample_id}.{library_id}.otu_table.tsv"
            for sample_id, library_id in SAMPLE_LIB
        ],


rule pre_singlem_condense:
    """Aggregate all the singlem results into a single table"""
    input:
        archive_otu_tables=[
            SINGLEM / f"{sample_id}.{library_id}.archive.json"
            for sample_id, library_id in SAMPLE_LIB
        ],
        data=rules.pre_singlem_data.output,
    output:
        condense=SINGLEM / "singlem.tsv",
    log:
        SINGLEM / "singlem.log",
    conda:
        "../envs/pre.yml"
    params:
        input_dir=SINGLEM,
    shell:
        """
        singlem condense \
            --input-archive-otu-tables {input.archive_otu_tables} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
        2> {log} 1>&2
        """


rule pre_singlem:
    """Run all stats singlem rules"""
    input:
        SINGLEM / "singlem.tsv",


rule pre_cram_to_mapped_bam:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=BOWTIE2 / "{sample}.{library_id}.cram",
        reference=features["host"]["fasta"],
    output:
        bam=temp(COVERM / "bams/{sample}.{library_id}.bam"),
    log:
        COVERM / "bams/{sample}.{library_id}.log",
    conda:
        "../envs/pre.yml"
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


rule pre_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=COVERM / "bams/{sample}.{library_id}.bam",
    output:
        tsv=touch(COVERM / "genome/{sample}.{library_id}.tsv"),
    conda:
        "../envs/pre.yml"
    log:
        COVERM / "genome/{sample}.{library_id}.log",
    params:
        methods=params["coverm"]["genome"]["methods"],
        min_covered_fraction=params["coverm"]["genome"]["min_covered_fraction"],
        separator=params["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output} 2> {log} || true \
        """


rule pre_coverm:
    """Aggregate all the nonpareil results into a single table"""
    input:
        [
            COVERM / f"genome/{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIB
        ],
    output:
        COVERM / "coverm.tsv",
    log:
        COVERM / "coverm.log",
    conda:
        "../envs/pre.yml"
    params:
        input_dir=COVERM / "genome",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule pre:
    input:
        rules.pre_fastp_trim_all.input,
        rules.pre_fastp_fastqc_all.input,
        rules.pre_bowtie2_map_host_all.input,
        rules.pre_bowtie2_extract_nonhost_all.input,
        rules.pre_nonpareil.output,
        rules.pre_singlem.input,
        rules.pre_coverm.output,