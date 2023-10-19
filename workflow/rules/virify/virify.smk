# include: "functions.smk"
include: "run.smk"


# include: "eval.smk"


rule virify:
    input:
        rules.virify_all.input,
