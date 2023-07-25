rule assembly_megahit_one:
    input:
        forwards=get_forwards_for_megahit,
        reverses=get_reverses_for_megahit,
    output:
        assembly_folder=directory(MEGAHIT / "{sample_id}"),
        assembly_fasta=MEGAHIT / "{sample_id}.contigs.fa",
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
            -1 {params.forwards} \
            -2 {params.reverses} \
            {params.extra} \
            --out-dir {output.assembly_folder} \
        2> {log} 1>&2

        cp \
            {output.assembly_folder}/final.contigs.fa \
            {output.assembly_fasta} \
        2>> {log} 1>&2
        """


rule assembly_megahit_all:
    input:
        [MEGAHIT / f"{sample_id}.contigs.fa" for sample_id in SAMPLES],


rule assembly:
    input:
        rules.assembly_megahit_all.input,
