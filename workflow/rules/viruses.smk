include: "viruses/mvp.smk"
include: "viruses/annotate.smk"
include: "viruses/multiqc.smk"


rule viruses__all:
    input:
        rules.viruses__mvp__all.input,
        rules.viruses__annotate__all.input,
        rules.viruses__multiqc__all.input,
