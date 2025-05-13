use rule bowtie2__build as prokaryotes__quantify__bowtie2__build with:
    input:
        ref=PROK_ANN / "drep.{secondary_ani}.fa.gz",
    output:
        multiext(
            str(PROK_BUILD / "drep.{secondary_ani}."),
            "1.bt2",
            "2.bt2",
            "3.bt2",
            "4.bt2",
            "rev.1.bt2",
            "rev.2.bt2",
        ),
    log:
        PROK_BUILD / "drep.{secondary_ani}.log",


rule prokaryotes__quantify__bowtie2__build__all:
    input:
        [
            PROK_BUILD / f"drep.{secondary_ani}.{extension}"
            for secondary_ani in SECONDARY_ANIS
            for extension in [
                "1.bt2",
                "2.bt2",
                "3.bt2",
                "4.bt2",
                "rev.1.bt2",
                "rev.2.bt2",
            ]
        ],


use rule bowtie2__map as prokaryotes__quantify__bowtie2__map with:
    input:
        forward_=PRE_CLEAN / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_CLEAN / "{sample_id}.{library_id}_2.fq.gz",
        mock=multiext(
            str(PROK_BUILD / "drep.{secondary_ani}."),
            "1.bt2",
            "2.bt2",
            "3.bt2",
            "4.bt2",
            "rev.1.bt2",
            "rev.2.bt2",
        ),
    output:
        PROK_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.bam",
    log:
        PROK_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.log",
    conda:
        "../../../environments/bowtie2.yml"
    params:
        index=lambda w: PROK_BUILD / f"drep.{w.secondary_ani}",
        bowtie2_extra=params["preprocess"]["bowtie2"]["bowtie2_extra"],
        samtools_extra=params["preprocess"]["bowtie2"]["samtools_extra"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,


rule prokaryotes__quantify__bowtie2__map__all:
    """Align all samples to the dereplicated genomes"""
    input:
        [
            PROK_BOWTIE2 / f"drep.{secondary_ani}" / f"{sample_id}.{library_id}.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
            for secondary_ani in SECONDARY_ANIS
        ],


rule prokaryotes__quantify__bowtie2__all:
    input:
        rules.prokaryotes__quantify__bowtie2__build__all.input,
        rules.prokaryotes__quantify__bowtie2__map__all.input,
