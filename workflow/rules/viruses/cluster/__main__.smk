include: "genomad.smk"
include: "dedupe.smk"
include: "mmseqs.smk"


rule viruses__cluster__all:
    input:
        rules.viruses__cluster__genomad__all.input,
        rules.viruses__cluster__dedupe__all.input,
        rules.viruses__cluster__mmseqs__all.input,
