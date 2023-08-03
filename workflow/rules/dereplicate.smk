include: "dereplicate_run.smk"
include: "dereplicate_eval.smk"


rule dereplicate:
    input:
        DREP,
        rules.dereplicate_eval_coverm_genome.output,
