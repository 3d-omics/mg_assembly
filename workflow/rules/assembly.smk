rule assembly_megahit_one:
    """Run megahit over one sample, merging all libraries in the process

    Note: the initial rm -rf is to delete the folder that snakemake creates.
    megahit refuses to overwrite an existing folder
    """
    input:
        forwards=get_forwards_for_megahit,
        reverses=get_reverses_for_megahit,
    output:
        assembly_folder=directory(MEGAHIT / "{sample_id}"),
        assembly_fasta=MEGAHIT / "{sample_id}/final.contigs.fa",
    log:
        log=MEGAHIT / "{sample_id}.log",
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
        [MEGAHIT / f"{sample_id}.contigs.fa" for sample_id in SAMPLES],


rule assembly_quast_one:
    """Run quast over one assembly group"""
    input:
        MEGAHIT / "{sample_id}/final.contigs.fa",
    output:
        directory(QUAST / "{sample_id}"),
    log:
        QUAST / "{sample_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 4
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """


rule assembly_quast_all:
    """Run quast over all assembly groups"""
    input:
        [QUAST / f"{sample_id}" for sample_id in SAMPLES],


rule assembly:
    """Run all the assemblies"""
    input:
        rules.assembly_quast_all.input,
