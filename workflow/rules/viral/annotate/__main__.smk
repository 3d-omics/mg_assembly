include: "dramv.smk"
include: "genomad.smk"
include: "quast.smk"
include: "virsorter2.smk"

    
rule viral__annotate:
    input:
        rules.viral__annotate__genomad.input,
        rules.viral__annotate__dramv.input,
        rules.viral__annotate__quast.input,
