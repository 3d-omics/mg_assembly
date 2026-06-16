include: "annotate/dramv.smk"
include: "annotate/genomad.smk"
include: "annotate/quast.smk"
include: "annotate/virsorter2.smk"
include: "annotate/checkv.smk"


rule viruses__annotate__all:
    input:
        rules.viruses__annotate__genomad__all.input,
        rules.viruses__annotate__dramv__all.input,
        rules.viruses__annotate__quast__all.input,
        rules.viruses__annotate__checkv__all.input,
