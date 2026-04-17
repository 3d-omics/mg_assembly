include: "megahit_functions.smk"


rule assemble__megahit:
    """Run megahit over one sample, merging all libraries in the process

Note: the initial rm -rf is to delete the folder that snakemake creates.
megahit refuses to overwrite an existing folder
"""
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        temp(directory(ASMB_MEGAHIT / "{assembly_id}.dir")),
    log:
        log=ASMB_MEGAHIT / "{assembly_id}.log",
    retries: 5
    conda:
        ENVS / "megahit.yml"
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=7 * 24 * 60,
    params:
        forwards=aggregate_forwards_for_megahit,
        reverses=aggregate_reverses_for_megahit,
        extra=params["assemble"]["megahit"]["extra"],
    shell:
        """
        megahit \
            --num-cpu-threads {threads} \
            --verbose \
            --force \
            --out-dir {output} \
            --continue \
            {params.extra} \
            -1 {params.forwards} \
            -2 {params.reverses} \
        2> {log} 1>&2
        """


rule assemble__megahit__rename:
    input:
        ASMB_MEGAHIT / "{assembly_id}.dir",
    output:
        ASMB_MEGAHIT / "{assembly_id}.fa.gz",
    log:
        ASMB_MEGAHIT / "{assembly_id}.rename.log",
    conda:
        "../../environments/megahit.yml"
    threads: 24
    params:
        assembly_id=lambda w: w.assembly_id,
    shell:
        """
        ( seqtk seq \
            {input}/final.contigs.fa \
        | cut -f 1 -d " " \
        | paste - - \
        | awk \
            '{{printf(">{params.assembly_id}:bin_NA@contig_%08d\\n%s\\n", NR, $2)}}' \
        | bgzip \
            --threads {threads} \
        > {output} \
        ) 2> {log}
        """


rule assemble__megahit__all:
    """Rename all assemblies contigs to avoid future collisions"""
    input:
        [ASMB_MEGAHIT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
