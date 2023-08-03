rule assemble_megahit_one:
    """Run megahit over one sample, merging all libraries in the process

    Note: the initial rm -rf is to delete the folder that snakemake creates.
    megahit refuses to overwrite an existing folder
    """
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        assemble_folder=directory(MEGAHIT / "{assemble_id}"),
        assemble_fasta=MEGAHIT / "{assemble_id}/final.contigs.fa",
    log:
        log=MEGAHIT / "{assemble_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 16
    resources:
        mem_mb=64 * 1024,
    params:
        min_contig_len=params["assemble"]["megahit"]["min_contig_len"],
        extra=params["assemble"]["megahit"]["extra"],
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
            --out-dir {output.assemble_folder} \
        2> {log} 1>&2
        """


rule assemble_megahit_all:
    """Run megahit over all groups"""
    input:
        [MEGAHIT / f"{assemble_id}/final.contigs.fa" for assemble_id in ASSEMBLIES],


rule assemble_megahit_renaming_one:
    """Rename contigs to avoid future collisions

    Contigs are renamed to `megahit_{assemble_id}.{contig_number}`.
    Also, secondary info in the header is removed.
    """
    input:
        MEGAHIT / "{assemble_id}/final.contigs.fa",
    output:
        ASSEMBLY_RENAME / "{assemble_id}.fa",
    log:
        ASSEMBLY_RENAME / "{assemble_id}.log",
    conda:
        "../envs/assembly.yml"
    params:
        assemble_id=lambda wildcards: wildcards.assemble_id,
    shell:
        """
        seqtk rename \
            {input} \
            {params.assemble_id}.contig \
        | cut -f 1 -d " " \
        > {output} 2> {log}
        """


rule assemble_bowtie2_build_one:
    """
    Index megahit assembly
    """
    input:
        contigs=ASSEMBLY_RENAME / "{assemble_id}.fa",
    output:
        mock=touch(ASSEMBLY_INDEX / "{assemble_id}"),
    log:
        ASSEMBLY_INDEX / "{assemble_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 24
    params:
        extra=params["assemble"]["bowtie2-build"]["extra"],
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule assemble_bowtie2_build_all:
    input:
        [ASSEMBLY_INDEX / f"{assemble_id}" for assemble_id in ASSEMBLIES],


rule assemble_bowtie2_one:
    input:
        mock=ASSEMBLY_INDEX / "{assemble_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        reference=ASSEMBLY_RENAME / "{assemble_id}.fa",
    output:
        cram=ASSEMBLY_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.cram",
    log:
        log=ASSEMBLY_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 24
    params:
        extra=params["assemble"]["bowtie2"]["extra"],
        samtools_mem=params["assemble"]["samtools"]["mem"],
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


rule assemble_bowtie2_all:
    input:
        [
            ASSEMBLY_BOWTIE2 / f"{assemble_id}.{sample_id}.{library_id}.cram"
            for assemble_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule assemble_cram_to_bam_one:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=ASSEMBLY_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.cram",
        crai=ASSEMBLY_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.cram.crai",
        reference=ASSEMBLY_RENAME / "{assemble_id}.fa",
    output:
        bam=temp(ASSEMBLY_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.bam"),
    log:
        ASSEMBLY_BOWTIE2 / "{assemble_id}.{sample_id}.{library_id}.bam.log",
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


rule assemble_cram_to_bam_all:
    input:
        [
            ASSEMBLY_BOWTIE2 / f"{assemble_id}.{sample_id}.{library_id}.bam"
            for assemble_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule assemble_run:
    input:
        rules.assemble_cram_to_bam_all.input,
