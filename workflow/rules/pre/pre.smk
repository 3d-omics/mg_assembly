include: "functions.smk"
include: "run.smk"
include: "eval.smk"


rule pre:
    input:
        rules.pre_run.input,
        rules.pre_eval.input,
