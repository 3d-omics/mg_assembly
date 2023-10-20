rule dereplicate_bowtie2_build_one:
    """
    Index dereplicader
    """
    input:
        contigs=DREP / "dereplicated_genomes.fa",
    output:
        mock=touch(DREP_INDEX / "dereplicated_genomes"),
    log:
        DREP_INDEX / "dereplicated_genomes.log",
    conda:
        "dereplicate.yml"
    threads: 24
    params:
        extra=params["dereplicate"]["bowtie2-build"]["extra"],
    resources:
        mem_mb=double_ram(params["dereplicate"]["bowtie2-build"]["memory_gb"]),
        runtime=24 * 60,
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule dereplicate_bowtie2_one:
    input:
        mock=DREP_INDEX / "dereplicated_genomes",
        forward_=get_forward_for_dereplicate_bowtie2_one,
        reverse_=get_reverse_for_dereplicate_bowtie2_one,
        reference=DREP / "dereplicated_genomes.fa",
    output:
        cram=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram",
    log:
        DREP_BOWTIE2 / "{sample_id}.{library_id}.log",
    conda:
        "dereplicate.yml"
    threads: 24
    params:
        extra=params["dereplicate"]["bowtie2"]["extra"],
        samtools_mem=params["dereplicate"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    resources:
        mem_mb=double_ram(params["dereplicate"]["bowtie2"]["memory_gb"]),
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


rule dereplicate_bowtie2:
    input:
        [
            DREP_BOWTIE2 / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],