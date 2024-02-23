include: "dramv.smk"
include: "genomad.smk"
include: "quast.smk"
include: "virsorter2.smk"


rule viruses__annotate:
    input:
        rules.viruses__annotate__genomad.output,
        rules.viruses__annotate__dramv.output,
        rules.viruses__annotate__quast.output,
