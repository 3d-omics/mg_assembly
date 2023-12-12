include: "__functions__.smk"


rule _reads__link:
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
        "__env__.yml"
    resources:  # run it superfast
        mem_mb=10,
        runtime=1,
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_} 2>  {log} 1>&2
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_} 2>> {log} 1>&2
        """


rule reads__link:
    """Link all reads in the samples.tsv"""
    input:
        [
            READS / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule reads__fastqc:
    """Get all fastqc reports of the raw reads"""
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.zip"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule reads__eval:
    """Perform the evaluation steps of the reads"""
    input:
        rules.reads__fastqc.input,


rule reads__run:
    """Only perform the "run" in reads without the evaluation steps"""
    input:
        rules.reads__link.input,


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads__run.input,
        rules.reads__eval.input,
