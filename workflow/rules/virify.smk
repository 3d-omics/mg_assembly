# include: "virify/functions.smk"
include: "virify/run.smk"


# include: "virify/eval.smk"


rule virify:
    input:
        rules.virify_all.input,
