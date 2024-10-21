include: "cluster/__main__.smk"
include: "annotate/__main__.smk"
include: "quantify/__main__.smk"


rule viruses__all:
    input:
        rules.viruses__cluster__all.input,
        rules.viruses__annotate__all.input,
        rules.viruses__quantify__all.input,
