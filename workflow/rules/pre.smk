include: "pre_run.smk"
include: "pre_eval.smk"


rule pre:
    input:
        rules.pre_run.input,
        rules.pre_eval.input,
