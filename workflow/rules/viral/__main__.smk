include: "__functions__.smk"
include: "genomad.smk"
include: "checkv.smk"
include: "diamond.smk"


rule viral:
    input:
        rules.viral__checkv.input,
        rules.viral__diamond.input,
