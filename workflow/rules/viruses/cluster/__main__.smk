include: "__functions__.smk"
include: "genomad.smk"
include: "dedupe.smk"
include: "mmseqs.smk"


rule viruses__cluster:
    input:
        rules.viruses__cluster__genomad.input,
        rules.viruses__cluster__dedupe.input,
        rules.viruses__cluster__mmseqs.input,
