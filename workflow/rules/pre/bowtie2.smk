rule pre_bowtie2_hosts_build:
    """Build PRE_BOWTIE2 index for the host reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "{genome}.fa.gz",
        faidx=REFERENCE / "{genome}.fa.gz.fai",
    output:
        mock=touch(PRE_INDEX / "{genome}"),
    log:
        PRE_INDEX / "{genome}.log",
    conda:
        "pre.yml"
    params:
        extra=params["pre"]["bowtie2-build"]["extra"],
    threads: 24
    resources:
        mem_mb=double_ram(params["pre"]["bowtie2-build"]["memory_gb"]),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.reference} \
            {output.mock} \
        2> {log} 1>&2
        """


rule pre_bowtie2_host_map_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_input_forward_for_host_mapping,
        reverse_=get_input_reverse_for_host_mapping,
        mock=PRE_INDEX / "{genome}",
        reference=REFERENCE / "{genome}.fa.gz",
        faidx=REFERENCE / "{genome}.fa.gz.fai",
    output:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
    log:
        PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.log",
    params:
        extra=params["pre"]["bowtie2"]["extra"],
        samtools_mem=params["pre"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    conda:
        "pre.yml"
    resources:
        mem_mb=double_ram(params["pre"]["bowtie2"]["memory_gb"]),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        find \
            $(dirname {output.cram}) \
            -name "$(basename {output.cram}).tmp.*.bam" \
            -delete \
        2> {log} 1>&2

        ( bowtie2 \
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
        ) 2>> {log} 1>&2
        """


rule pre_bowtie2_extract_nonhost_one:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
        reference=REFERENCE / "{genome}.fa.gz",
        fai=REFERENCE / "{genome}.fa.gz.fai",
    output:
        forward_=PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}_2.fq.gz",
    log:
        PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}.log",
    conda:
        "pre.yml"
    params:
        samtools_mem=params["pre"]["bowtie2"]["samtools"]["mem_per_thread"],
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=32 * 1024,
    shell:
        """
        ( samtools view \
            --reference {input.reference} \
            --threads {threads} \
            -u \
            -o /dev/stdout \
            -f 12 \
            {input.cram} \
        | samtools collate \
            -O \
            -u \
            -f \
            --reference {input.reference} \
            -@ {threads} \
            - \
        | samtools fastq \
            -1 >(pigz > {output.forward_}) \
            -2 >(pigz > {output.reverse_}) \
            -0 /dev/null \
            -c 9 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule pre_bowtie2_extract_nonhost_all:
    """Run bowtie2_extract_nonchicken_one for all PE libraries"""
    input:
        [
            PRE_BOWTIE2 / f"non{genome}" / f"{sample_id}.{library_id}_{end}.fq.gz"
            for genome in [LAST_HOST]
            if LAST_HOST
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule pre_bowtie2:
    """Run all the preprocessing steps for bowtie2"""
    input:
        rules.pre_bowtie2_extract_nonhost_all.input,
