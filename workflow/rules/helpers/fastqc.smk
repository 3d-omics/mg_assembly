rule fastqc:
    input:
        "{prefix}.fq.gz",
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip",
    conda:
        "helpers.yml"
    log:
        "{prefix}_fastqc.log",
    shell:
        "fastqc {input} 2> {log} 1>&2"
