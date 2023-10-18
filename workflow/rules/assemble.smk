include: "assemble/functions.smk"
include: "assemble/run.smk"
include: "assemble/eval.smk"


rule assemble:
    """Run the assemble module"""
    input:
        rules.assemble_eval.input,
        rules.assemble_run.input,
