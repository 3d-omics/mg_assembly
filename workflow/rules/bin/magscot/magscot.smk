include: "functions.smk"
include: "run.smk"


rule magscot:
    input:
        rules.metabin_run.input,
