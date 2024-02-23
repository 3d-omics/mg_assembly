include: "cluster/__main__.smk"
include: "annotate/__main__.smk"
include: "quantify/__main__.smk"

rule prokaryotes:
    input:
        rules.prokaryotes__cluster.input,
        rules.prokaryotes__annotate.input,
        rules.prokaryotes__quantify.input,
