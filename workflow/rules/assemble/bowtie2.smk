use rule bowtie2__build as assemble__bowtie2__build with:
    input:
        ref=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
    output:
        multiext(
            str(ASMB_BUILD / "{assembly_id}"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        ASMB_BUILD / "{assembly_id}.log",
    retries: 5


rule assemble__bowtie2__build__all:
    """Index all megahit assemblies"""
    input:
        [
            ASMB_BUILD / f"{assembly_id}.{extension}"
            for assembly_id in ASSEMBLIES
            for extension in [
                "1.bt2",
                "2.bt2",
                "3.bt2",
                "4.bt2",
                "rev.1.bt2",
                "rev.2.bt2",
            ]
        ],


use rule bowtie2__map as assemble__bowtie2__map with:
    input:
        forward_=PRE_CLEAN / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_CLEAN / "{sample_id}.{library_id}_2.fq.gz",
        mock=multiext(
            str(ASMB_BUILD / "{assembly_id}"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        ASMB_BOWTIE2 / "{assembly_id}" / "{sample_id}.{library_id}.bam",
    log:
        ASMB_BOWTIE2 / "{assembly_id}" / "{sample_id}.{library_id}.log",
    params:
        index=lambda w: ASMB_BUILD / f"{w.assembly_id}",
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        samtools_extra="",
        bowtie2_extra="",


rule assemble__bowtie2__map__all:
    """Map all samples to all the assemblies that they belong to"""
    input:
        [
            ASMB_BOWTIE2 / assembly_id / f"{sample_id}.{library_id}.bam"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule assemble__bowtie2__all:
    input:
        rules.assemble__bowtie2__build__all.input,
        rules.assemble__bowtie2__map__all.input,
