include: "reads_functions.smk"
include: "reads_run.smk"
include: "reads_eval.smk"


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads_run.input,
        rules.reads_eval.input,
