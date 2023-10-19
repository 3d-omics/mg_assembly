include: "functions.smk"
include: "run.smk"
include: "eval.smk"


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads_run.input,
        rules.reads_eval.input,
