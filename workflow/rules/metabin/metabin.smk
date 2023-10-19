include: "functions.smk"
include: "run.smk"


rule metabin:
    input:
        rules.metabin_run.input,
