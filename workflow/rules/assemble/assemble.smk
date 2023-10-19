include: "functions.smk"
include: "run.smk"
include: "eval.smk"


rule assemble:
    """Run the assemble module"""
    input:
        rules.assemble_eval.input,
        rules.assemble_run.input,
