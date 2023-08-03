include: "metabin/magscot.smk"


rule metabin_index_one:
    input:
        bins=MAGSCOT / "{assembly_id}.fa",
    output:
        mock=touch(METABIN_INDEX / "{assembly_id}"),
    log:
        METABIN_INDEX / "{assembly_id}.log",
    conda:
        "../envs/metabin.yml"
    threads: 24
    params:
        extra=params["metabin"]["bowtie2-build"]["extra"],
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.bins} \
            {output.mock} \
        2> {log} 1>&2
        """


rule metabin_bowtie2_one:
    input:
        mock=METABIN_INDEX / "{assembly_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        reference=MAGSCOT / "{assembly_id}.fa",
    output:
        cram=METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
    log:
        METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        "../envs/metabin.yml"
    threads: 24
    params:
        extra=params["metabin"]["bowtie2"]["extra"],
        samtools_mem=params["metabin"]["samtools"]["mem"],
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


rule metabin_bowtie2:
    input:
        [
            METABIN_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule metabin_run:
    input:
        rules.metabin_bowtie2.input,
