include: "__functions__.smk"
include: "genomad.smk"
include: "checkv.smk"
include: "diamond.smk"
include: "htseq.smk"
include: "virsorter.smk"
include: "dramv.smk"
include: "dedupe.smk"
include: "mmseqs.smk"


rule viral:
    input:
        rules.viral__checkv.input,
        rules.viral__diamond.input,
        rules.viral__htseq.input,
        rules.viral__virsorter.input,
        rules.viral__dramv.input,
        rules.viral__dedupe.input,
        rules.viral__mmseqs.input,
