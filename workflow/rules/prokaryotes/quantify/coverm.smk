use rule coverm__genome as prokaryotes__quantify__coverm__genome with:
    input:
        PROK_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.bam",
    output:
        temp(
            PROK_COVERM
            / "genome"
            / "{method}.drep.{secondary_ani}.{sample_id}.{library_id}.tsv.gz"
        ),
    log:
        PROK_COVERM
        / "genome"
        / "{method}.drep.{secondary_ani}.{sample_id}.{library_id}.log",
    params:
        method=lambda w: w.method,
        extra=params["quantify"]["coverm"]["genome"]["extra"],
        separator=params["quantify"]["coverm"]["genome"]["separator"],


use rule csvtk__join__left as prokaryotes__quantify__coverm__genome__join with:
    input:
        lambda w: [
            PROK_COVERM
            / "genome"
            / f"{w.method}.drep.{w.secondary_ani}.{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
        + ["/dev/null"],
    output:
        PROK_COVERM / "genome.{method}.drep.{secondary_ani}.tsv.gz",
    log:
        PROK_COVERM / "genome.{method}.drep.{secondary_ani}.log",


rule prokaryotes__quantify__coverm__genome__all:
    """Run coverm genome and all methods"""
    input:
        [
            PROK_COVERM / f"genome.{method}.drep.{secondary_ani}.tsv.gz"
            for method in ["count", "covered_bases"]
            for secondary_ani in SECONDARY_ANIS
        ],


rule prokaryotes__quantify__coverm__all:
    input:
        rules.prokaryotes__quantify__coverm__genome__all.input,
