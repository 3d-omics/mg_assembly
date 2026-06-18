rule viruses__mvp__04_read_mapping:
    """Map one sample's reads against the representative viral catalog, in parallel across samples"""
    input:
        metadata=VIR_MVP / "metadata.tsv",
        bowtie2_index=multiext(
            str(VIR_MVP_READ_MAPPING / "reference"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        bam=VIR_MVP_READ_MAPPING / "{sample_id}" / "{sample_id}_sorted.bam",
        coverm=VIR_MVP_READ_MAPPING / "{sample_id}" / "{sample_id}_CoverM.tsv",
    log:
        VIR_MVP_READ_MAPPING / "{sample_id}.log",
    conda:
        ENVS / "mvp.yml"
    threads: 24
    resources:
        mem_mb=double_ram(16 * 1024),
        runtime=6 * 60,
    params:
        sample_number=get_sample_number,
        read_type=params["viral"]["mvp"]["read_mapping"]["read_type"],
        interleaved=params["viral"]["mvp"]["read_mapping"]["interleaved"],
    shell:
        """
        mvip MVP_04_do_read_mapping \
            --input {VIR_MVP} \
            --metadata {input.metadata} \
            --sample_group {params.sample_number} \
            --read_type {params.read_type} \
            --interleaved {params.interleaved} \
            --threads {threads} \
        2> {log} 1>&2
        """


rule viruses__mvp__04_read_mapping__all:
    input:
        [
            VIR_MVP_READ_MAPPING / f"{sample_id}" / f"{sample_id}_CoverM.tsv"
            for sample_id in SAMPLES
        ],
