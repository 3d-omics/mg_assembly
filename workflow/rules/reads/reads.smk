include: "functions.smk"


rule reads_link_one:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=READS / "{sample}.{library}_1.fq.gz",
        reverse_=READS / "{sample}.{library}_2.fq.gz",
    log:
        READS / "{sample}.{library}.log",
    conda:
        "reads.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_} 2>  {log} 1>&2
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_} 2>> {log} 1>&2
        """


rule reads_link_all:
    """Link all reads in the samples.tsv"""
    input:
        [
            READS / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule reads_fastqc:
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.zip"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule reads_eval:
    input:
        rules.reads_fastqc.input,


rule reads_run:
    input:
        rules.reads_link_all.input,


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads_run.input,
        rules.reads_eval.input,


localrules:
    reads_link_one,
