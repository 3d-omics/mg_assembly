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
    conda:
        "../../../environments/coverm.yml"
    params:
        method=lambda w: w.method,
        extra=params["quantify"]["coverm"]["genome"]["extra"],
        separator=params["quantify"]["coverm"]["genome"]["separator"],


rule prokaryotes__quantify__coverm__genome__join:
    input:
        lambda w: [
            PROK_COVERM
            / "genome"
            / f"{w.method}.drep.{w.secondary_ani}.{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        PROK_COVERM / "genome.{method}.drep.{secondary_ani}.tsv.gz",
    log:
        PROK_COVERM / "genome.{method}.drep.{secondary_ani}.log",
    params:
        subcommand="join",
        extra="--left-join --tabs --out-tabs",
    resources:
        runtime=60,
        mem_mb=8 * 1024,
    wrapper:
        "v5.2.1/utils/csvtk"


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
