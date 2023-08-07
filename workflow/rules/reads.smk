include: "reads/functions.smk"
include: "reads/run.smk"
include: "reads/eval.smk"


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads_run.input,
        rules.reads_eval.input,
