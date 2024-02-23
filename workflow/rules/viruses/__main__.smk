include: "cluster/__main__.smk"
include: "annotate/__main__.smk"
include: "quantify/__main__.smk"


rule viruses:
    input:
        rules.viruses__cluster.input,
        rules.viruses__annotate.input,
        rules.viruses__quantify.input,
