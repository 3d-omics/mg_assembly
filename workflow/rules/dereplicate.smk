include: "dereplicate/run.smk"
include: "dereplicate/eval.smk"


rule dereplicate:
    input:
        rules.dereplicate_run.input,
        rules.dereplicate_eval.input,
