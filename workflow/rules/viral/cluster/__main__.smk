include: "__functions__.smk"
include: "genomad.smk"
include: "dedupe.smk"
include: "mmseqs.smk"


rule viral__cluster:
    input:
        rules.viral__cluster__genomad.input,
        rules.viral__cluster__dedupe.input,
        rules.viral__cluster__mmseqs.input,
