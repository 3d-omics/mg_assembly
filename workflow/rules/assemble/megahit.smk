rule assemble_megahit_one:
    """Run megahit over one sample, merging all libraries in the process

    Note: the initial rm -rf is to delete the folder that snakemake creates.
    megahit refuses to overwrite an existing folder
    """
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        fasta=MEGAHIT / "{assembly_id}.fa",
        tarball=MEGAHIT / "{assembly_id}.tar.gz",
    log:
        log=MEGAHIT / "{assembly_id}.log",
    conda:
        "assemble.yml"
    threads: 24
    params:
        out_dir=compose_out_dir_for_assemble_megahit_one,
        min_contig_len=params["assemble"]["megahit"]["min_contig_len"],
        extra=params["assemble"]["megahit"]["extra"],
        forwards=aggregate_forwards_for_megahit,
        reverses=aggregate_reverses_for_megahit,
    resources:
        mem_mb=double_ram(params["assemble"]["megahit"]["memory_gb"]),
        runtime=7 * 24 * 60,
    retries: 5
    shell:
        """
        megahit \
            --num-cpu-threads {threads} \
            --min-contig-len {params.min_contig_len} \
            --verbose \
            --force \
            --out-dir {params.out_dir} \
            --continue \
            -1 {params.forwards} \
            -2 {params.reverses} \
            {params.extra} \
        2> {log} 1>&2

        cp {params.out_dir}/final.contigs.fa {output.fasta} 2>> {log} 1>&2

        tar \
            --create \
            --file {output.tarball} \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            {params.out_dir} \
        2>> {log} 1>&2
        """


rule assemble_megahit_all:
    """Run megahit over all groups"""
    input:
        [MEGAHIT / f"{assembly_id}.fa" for assembly_id in ASSEMBLIES],


rule assemble_megahit_renaming_one:
    """Rename contigs to avoid future collisions

    Contigs are renamed to `megahit_{assembly_id}.{contig_number}`.
    Also, secondary info in the header is removed.
    """
    input:
        MEGAHIT / "{assembly_id}.fa",
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


rule assemble_megahit:
    """Rename all assemblies contigs to avoid future collisions"""
    input:
        [ASSEMBLE_RENAME / f"{assembly_id}.fa" for assembly_id in ASSEMBLIES],
