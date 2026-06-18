include: "annotate/dramv.smk"
include: "annotate/quast.smk"


rule viruses__annotate__all:
    input:
        rules.viruses__annotate__dramv__all.input,
        rules.viruses__annotate__quast__all.input,
