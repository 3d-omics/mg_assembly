# coverm contig ----
use rule coverm__contig as viruses__quantify__coverm__contig with:
    input:
        VIR_BOWTIE2 / "rep_seq" / "{sample_id}.{library_id}.bam",
    output:
        temp(VIR_COVERM / "contig" / "{method}.rep_seq.{sample_id}.{library_id}.tsv.gz"),
    log:
        VIR_COVERM / "contig" / "{method}.{sample_id}.{library_id}.log",
    params:
        method=lambda w: w.method,


use rule csvtk__join__left as viruses__quantify__coverm__contig__join with:
    input:
        lambda w: [
            VIR_COVERM
            / "contig"
            / f"{w.method}.rep_seq.{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
        + ["/dev/null"],
    output:
        VIR_COVERM / "contig.{method}.rep_seq.tsv.gz",
    log:
        VIR_COVERM / "contig.{method}.rep_seq.log",


rule viruses__quantify__coverm__contig__all:
    """Run coverm contig and all methods"""
    input:
        [
            VIR_COVERM / f"contig.{method}.rep_seq.tsv.gz"
            for method in ["count", "covered_bases"]
        ],


rule viruses__quantify__coverm__all:
    input:
        rules.viruses__quantify__coverm__contig__all.input,
