include: "__functions__.smk"
include: "genomad.smk"


rule viral:
    input:
        rules.viral__genomad.input,
