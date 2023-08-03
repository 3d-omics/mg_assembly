include: "assemble_run.smk"
include: "assemble_eval.smk"


rule assemble:
    input:
        rules.assemble_eval.input,
        rules.assemble_run.input,
