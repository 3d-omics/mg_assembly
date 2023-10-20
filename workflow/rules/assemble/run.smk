rule assemble_megahit_one:
    """Run megahit over one sample, merging all libraries in the process

    Note: the initial rm -rf is to delete the folder that snakemake creates.
    megahit refuses to overwrite an existing folder
    """
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        # assemble_folder=directory(MEGAHIT / "{assembly_id}"),
        assemble_fasta=MEGAHIT / "{assembly_id}/final.contigs.fa",
    log:
        log=MEGAHIT / "{assembly_id}.log",
    conda:
        "assemble.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["assemble"]["megahit"]["memory_gb"]),
    params:
        out_dir=compose_out_dir_for_assemble_megahit_one,
        min_contig_len=params["assemble"]["megahit"]["min_contig_len"],
        extra=params["assemble"]["megahit"]["extra"],
        forwards=aggregate_forwards_for_megahit,
        reverses=aggregate_reverses_for_megahit,
        memory_bytes=get_memory_bytes_for_assemble_megahit_one,
    retries: 5
    shell:
        """
        megahit \
            --num-cpu-threads {threads} \
            --min-contig-len {params.min_contig_len} \
            --memory {params.memory_bytes} \
            --verbose \
            --force \
            --out-dir {params.out_dir} \
            -1 {params.forwards} \
            -2 {params.reverses} \
            {params.extra} \
        2> {log} 1>&2
        """


rule assemble_megahit_all:
    """Run megahit over all groups"""
    input:
        [MEGAHIT / f"{assembly_id}/final.contigs.fa" for assembly_id in ASSEMBLIES],


rule assemble_megahit_renaming_one:
    """Rename contigs to avoid future collisions

    Contigs are renamed to `megahit_{assembly_id}.{contig_number}`.
    Also, secondary info in the header is removed.
    """
    input:
        MEGAHIT / "{assembly_id}/final.contigs.fa",
    output:
        ASSEMBLE_RENAME / "{assembly_id}.fa",
    log:
        ASSEMBLE_RENAME / "{assembly_id}.log",
    conda:
        "assemble.yml"
    params:
        assembly_id="{assembly_id}",
    shell:
        """
        ( seqtk seq {input} \
        | cut -f 1 -d " " \
        | paste - - \
        | awk '{{printf(">{params.assembly_id}:bin_NA@contig_%08d\\n%s\\n", NR, $2)}}' \
        > {output} \
        ) 2> {log}
        """


rule assemble_bowtie2_build_one:
    """Index a megahit assembly"""
    input:
        contigs=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        mock=touch(ASSEMBLE_INDEX / "{assembly_id}"),
    log:
        ASSEMBLE_INDEX / "{assembly_id}.log",
    conda:
        "assemble.yml"
    threads: 24
    params:
        extra=params["assemble"]["bowtie2-build"]["extra"],
    resources:
        mem_mb=double_ram(params["assemble"]["bowtie2-build"]["memory_gb"]),
        runtime=6 * 60,
    retries: 5
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
    """Index all megahit assemblies"""
    input:
        [ASSEMBLE_INDEX / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule assemble_bowtie2_one:
    """Map one sample to one megahit assembly"""
    input:
        mock=ASSEMBLE_INDEX / "{assembly_id}",
        forward_=get_forwards_from_assembly_id,
        reverse_=get_reverses_from_assembly_id,
        reference=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        cram=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
    log:
        log=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        "assemble.yml"
    threads: 24
    params:
        extra=params["assemble"]["bowtie2"]["extra"],
        samtools_mem=params["assemble"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    resources:
        mem_mb=double_ram(params["assemble"]["bowtie2"]["memory_gb"]),
        runtime=6 * 60,
    retries: 5
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
    """Map all samples to all the assemblies that they belong to"""
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule assemble_run:
    """Run all assembly rules (no evaluation)"""
    input:
        rules.assemble_bowtie2_all.input,
