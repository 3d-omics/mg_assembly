include: "checkv.smk"
include: "virsorter2.smk"
include: "dramv.smk"
include: "quast.smk"


rule viral__annotate:
    input:
        rules.viral__annotate__checkv.input,
        rules.viral__annotate__virsorter2.input,
        rules.viral__annotate__dramv.input,
        rules.viral__annotate__quast.input,
