include: "functions.smk"
include: "run.smk"
include: "eval.smk"


rule dereplicate:
    input:
        rules.dereplicate_run.input,
        rules.dereplicate_eval.input,
