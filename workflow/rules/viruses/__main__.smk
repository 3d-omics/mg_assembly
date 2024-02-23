include: "cluster/__main__.smk"
include: "annotate/__main__.smk"
include: "quantify/__main__.smk"


rule viral:
    input:
        rules.viral__cluster.input,
        rules.viral__annotate.input,
        rules.viral__quantify.input,
