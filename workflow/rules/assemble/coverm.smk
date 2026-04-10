use rule coverm__contig as assemble__coverm__contig with:
    input:
        ASMB_BOWTIE2 / "{assembly_id}" / "{sample_id}.{library_id}.bam",
    output:
        temp(
            ASMB_COVERM
            / "contig"
            / "{method}.{assembly_id}.{sample_id}.{library_id}.tsv.gz"
        ),
    log:
        ASMB_COVERM / "contig" / "{method}.{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        ENVS / "coverm.yml"
    params:
        method=lambda w: w.method,


use rule csvtk__join__left as assemble__coverm__join with:
    input:
        lambda w: [
            ASMB_COVERM
            / "contig"
            / f"{w.method}.{w.assembly_id}.{sample_id}.{library_id}.tsv.gz"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
            if assembly_id == w.assembly_id
        ]
        + NULL,
    output:
        ASMB_COVERM / "contig.{method}.{assembly_id}.tsv.gz",
    log:
        ASMB_COVERM / "contig.{method}.{assembly_id}.log",


rule assemble__coverm__all:
    input:
        [
            ASMB_COVERM / f"contig.{method}.{assembly_id}.tsv.gz"
            for assembly_id in ASSEMBLIES
            for method in ["count", "covered_bases"]
        ],
