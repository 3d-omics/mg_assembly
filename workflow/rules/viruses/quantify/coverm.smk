# coverm contig ----
use rule coverm__contig as viruses__quantify__coverm__contig with:
    input:
        VIR_BOWTIE2 / "rep_seq" / "{sample_id}.{library_id}.bam",
    output:
        temp(VIR_COVERM / "contig" / "{method}.rep_seq.{sample_id}.{library_id}.tsv.gz"),
    log:
        VIR_COVERM / "contig" / "{method}.{sample_id}.{library_id}.log",
    conda:
        "../../../environments/coverm.yml"
    params:
        method=lambda w: w.method,


rule viruses__quantify__coverm__contig__join:
    input:
        lambda w: [
            VIR_COVERM
            / "contig"
            / f"{w.method}.rep_seq.{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        VIR_COVERM / "contig.{method}.rep_seq.tsv.gz",
    log:
        VIR_COVERM / "contig.{method}.rep_seq.log",
    params:
        subcommand="join",
        extra="--left-join --tabs --out-tabs",
    wrapper:
        "v5.2.1/utils/csvtk"


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
