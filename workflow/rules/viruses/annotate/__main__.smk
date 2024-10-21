include: "dramv.smk"
include: "genomad.smk"
include: "quast.smk"
include: "virsorter2.smk"
include: "checkv.smk"


rule viruses__annotate__all:
    input:
        rules.viruses__annotate__genomad__all.input,
        rules.viruses__annotate__dramv__all.input,
        rules.viruses__annotate__quast__all.input,
        rules.viruses__annotate__checkv__all.input,
        rules.viruses__annotate__quast__all.input,
