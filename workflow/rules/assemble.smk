include: "assemble/functions.smk"
include: "assemble/run.smk"
include: "assemble/eval.smk"


rule assemble:
    input:
        rules.assemble_eval.input,
        rules.assemble_run.input,
