include: "__functions__.smk"
include: "genomad.smk"
include: "checkv.smk"


rule viral:
    input:
        rules.viral__genomad.input,
