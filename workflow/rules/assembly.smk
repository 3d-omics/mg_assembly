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
        [MEGAHIT / f"{assembly_id}.contigs.fa" for assembly_id in samples.assembly_id],


rule assembly_quast_one:
    """Run quast over one assembly group"""
    input:
        MEGAHIT / "{assembly_id}/final.contigs.fa",
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


rule assembly:
    """Run all the assemblies"""
    input:
        rules.assembly_quast_all.input,


rule assembly_bowtie2_build_one:
    """
    Index megahit assembly
    """
    input:
        contigs=MEGAHIT / "{assembly_id}" / "final.contigs.fa",
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
        reference=MEGAHIT / "{assembly_id}" / "final.contigs.fa",
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
