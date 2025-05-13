use rule bowtie2__build as viruses__quantify__bowtie2__build with:
    input:
        ref=VIR_MMSEQS / "rep_seq.fa.gz",
    output:
        multiext(
            str(VIR_BUILD / "viruses"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        VIR_BUILD / "virues.log",


rule viruses__quantify__bowtie2__build__all:
    input:
        rules.viruses__quantify__bowtie2__build.output,


use rule bowtie2__map as viruses__quantify__bowtie2__map with:
    input:
        forward_=PRE_CLEAN / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_CLEAN / "{sample_id}.{library_id}_2.fq.gz",
        mock=multiext(
            str(VIR_BUILD / "viruses"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        VIR_BOWTIE2 / "rep_seq" / "{sample_id}.{library_id}.bam",
    log:
        VIR_BOWTIE2 / "rep_seq" / "{sample_id}.{library_id}.log",
    params:
        index=VIR_BUILD / "viruses",
        samtools_extra=params["preprocess"]["bowtie2"]["samtools_extra"],
        bowtie2_extra=params["preprocess"]["bowtie2"]["bowtie2_extra"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    conda:
        "../../../environments/bowtie2.yml"


rule viruses__quantify__bowtie2__map__all:
    """Align all samples to the dereplicated genomes"""
    input:
        [
            VIR_BOWTIE2 / "rep_seq" / f"{sample_id}.{library_id}.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule viruses__quantify__bowtie2__all:
    input:
        rules.viruses__quantify__bowtie2__build__all.input,
        rules.viruses__quantify__bowtie2__map__all.input,
