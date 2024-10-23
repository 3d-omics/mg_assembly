include: "genomad.smk"
include: "bbmap.smk"
include: "mmseqs.smk"


rule viruses__cluster__all:
    input:
        rules.viruses__cluster__genomad__all.input,
        rules.viruses__cluster__bbmap__all.input,
        rules.viruses__cluster__mmseqs__all.input,
