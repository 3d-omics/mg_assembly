include: "dramv.smk"
include: "genomad.smk"
include: "quast.smk"
include: "virsorter2.smk"
include: "checkv.smk"


rule viruses__annotate:
    input:
        rules.viruses__annotate__genomad.input,
        rules.viruses__annotate__dramv.input,
        rules.viruses__annotate__quast.input,
        rules.viruses__annotate__checkv.input,
