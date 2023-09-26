include: "pre/functions.smk"
include: "pre/run.smk"
include: "pre/eval.smk"


rule pre:
    input:
        rules.pre_run.input,
        rules.pre_eval.input,
