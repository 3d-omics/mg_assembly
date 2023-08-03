include: "dereplicate_run.smk"
include: "dereplicate_eval.smk"


rule dereplicate:
    input:
        rules.dereplicate_run.input,
        rules.dereplicate_eval.input,
