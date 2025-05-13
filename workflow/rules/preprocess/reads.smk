include: "reads_functions.smk"


rule preprocess__reads:
    """Make a link to the original forward file, with a prettier name than default"""
    input:
        forward_=get_forward_filename,
        reverse_=get_reverse_filename,
    output:
        forward_=PRE_READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_READS / "{sample_id}.{library_id}_2.fq.gz",
    log:
        PRE_READS / "{sample_id}.{library_id}.log",
    conda:
        "base"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_} 2> {log}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_} 2> {log}
        """


rule preprocess__reads__all:
    """Link all reads in the samples.tsv"""
    input:
        [
            PRE_READS / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],
