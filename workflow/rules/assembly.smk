rule assembly_megahit_one:
    """Run megahit over one sample, merging all libraries in the process

    Note: the initial rm -rf is to delete the folder that snakemake creates.
    megahit refuses to overwrite an existing folder
    """
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        assembly_folder=directory(MEGAHIT / "{assembly_id}"),
        assembly_fasta=MEGAHIT / "{assembly_id}/final.contigs.fa",
    log:
        log=MEGAHIT / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 16
    resources:
        mem_mb=64 * 1024,
    params:
        min_contig_len=params["assembly"]["megahit"]["min_contig_len"],
        extra=params["assembly"]["megahit"]["extra"],
        forwards=aggregate_forwards_for_megahit,
        reverses=aggregate_reverses_for_megahit,
    shell:
        """
        megahit \
            --num-cpu-threads {threads} \
            --min-contig-len {params.min_contig_len} \
            --verbose \
            --force \
            -1 {params.forwards} \
            -2 {params.reverses} \
            {params.extra} \
            --out-dir {output.assembly_folder} \
        2> {log} 1>&2
        """


rule assembly_megahit_all:
    """Run megahit over all groups"""
    input:
        [
            MEGAHIT / f"{assembly_id}/final.contigs.fa"
            for assembly_id in samples.assembly_id
        ],


rule assembly_megahit_rename_one:
    """Rename contigs to avoid future collisiions to `megahit_{assembly_id}.{contig_number}`"""
    input:
        MEGAHIT / "{assembly_id}/final.contigs.fa",
    output:
        MEGAHIT_RENAMED / "{assembly_id}.fa",
    log:
        MEGAHIT_RENAMED / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
    params:
        assembly_id=lambda wildcards: wildcards.assembly_id,
    shell:
        """
        seqtk rename \
            {input} \
            megahit_{params.assembly_id}. \
        > {output} 2> {log}
        """


rule assembly_quast_one:
    """Run quast over one assembly group"""
    input:
        MEGAHIT_RENAMED / "{assembly_id}.fa",
    output:
        directory(QUAST / "{assembly_id}"),
    log:
        QUAST / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 4
    params:
        extra=params["assembly"]["quast"]["extra"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """


rule assembly_quast_all:
    """Run quast over all assembly groups"""
    input:
        [QUAST / f"{assembly_id}" for assembly_id in samples.assembly_id],


rule assembly_bowtie2_build_one:
    """
    Index megahit assembly
    """
    input:
        contigs=MEGAHIT_RENAMED / "{assembly_id}.fa",
    output:
        mock=touch(BOWTIE2_INDEXES_ASSEMBLY / "{assembly_id}"),
    log:
        BOWTIE2_INDEXES_ASSEMBLY / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 24
    params:
        extra=params["assembly"]["bowtie2-build"]["extra"],
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule assembly_bowtie2_build_all:
    input:
        [
            BOWTIE2_INDEXES_ASSEMBLY / f"{assembly_id}"
            for assembly_id in samples.assembly_id
        ],


rule assembly_bowtie2_one:
    input:
        mock=BOWTIE2_INDEXES_ASSEMBLY / "{assembly_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        reference=MEGAHIT_RENAMED / "{assembly_id}.fa",
    output:
        cram=BOWTIE2_ASSEMBLY / "{assembly_id}/{sample_id}.{library_id}.cram",
    log:
        log=BOWTIE2_ASSEMBLY / "{assembly_id}/{sample_id}.{library_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 24
    params:
        extra=params["assembly"]["bowtie2"]["extra"],
        samtools_mem=params["assembly"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    resources:
        mem_mb=32 * 1024,
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


rule assembly_bowtie2_all:
    input:
        [
            BOWTIE2_ASSEMBLY / f"{assembly_id}/{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in samples[
                ["assembly_id", "sample_id", "library_id"]
            ].values.tolist()
        ],


rule assembly_cram_to_mapped_bam:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=BOWTIE2_ASSEMBLY / "{assembly_id}/{sample_id}.{library_id}.cram",
        reference=MEGAHIT_RENAMED / "{assembly_id}.fa",
    output:
        bam=temp(COVERM_ASSEMBLY / "bams/{assembly_id}/{sample_id}.{library_id}.bam"),
    log:
        COVERM_ASSEMBLY / "bams/{assembly_id}/{sample_id}.{library_id}.log",
    conda:
        "../envs/assembly.yml"
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


rule assembly_coverm_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=COVERM_ASSEMBLY / "bams/{assembly_id}/{sample_id}.{library_id}.bam",
        reference=MEGAHIT_RENAMED / "{assembly_id}.fa",
    output:
        tsv=COVERM_ASSEMBLY / "contig/{assembly_id}/{sample_id}.{library_id}.tsv",
    conda:
        "../envs/assembly.yml"
    log:
        COVERM_ASSEMBLY / "contig/{assembly_id}/{sample_id}.{library_id}.log",
    params:
        methods=params["assembly"]["coverm"]["genome"]["methods"],
        min_covered_fraction=params["assembly"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["assembly"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output.tsv} 2> {log}
        """


rule assembly_coverm_aggregate_one:
    input:
        tsvs=get_tsvs_for_assembly_coverm_genome,
    output:
        tsv=COVERM_ASSEMBLY / "{assembly_id}_contig.tsv",
    log:
        COVERM_ASSEMBLY / "{assembly_id}_contig.log",
    conda:
        "../envs/assembly.yml"
    params:
        input_dir=lambda wildcards: COVERM_ASSEMBLY / f"contig/{wildcards.assembly_id}",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule assembly_coverm_aggregate_all:
    input:
        tsvs=[
            COVERM_ASSEMBLY / f"{assembly_id}_contig.tsv"
            for assembly_id in samples.assembly_id
        ],


rule assembly:
    input:
        rules.assembly_coverm_aggregate_all.input,
