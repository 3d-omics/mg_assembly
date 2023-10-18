include: "metabin/functions.smk"
include: "metabin/run.smk"


rule metabin:
    input:
        rules.metabin_run.input,
