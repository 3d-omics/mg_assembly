# include: "metabin_functions.smk"
include: "metabin/run.smk"
include: "metabin/eval.smk"


rule metabin:
    input:
        rules.metabin_run.input,
        rules.metabin_eval.input,
