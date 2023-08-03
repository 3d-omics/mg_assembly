# include: "metabin_functions.smk"
include: "metabin_run.smk"
include: "metabin_eval.smk"


rule metabin:
    input:
        rules.metabin_run.input,
        rules.metabin_eval.input,
