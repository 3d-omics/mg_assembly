include: "mags.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "checkm2.smk"
include: "drep.smk"


rule prokaryotes__annotate:
    """Evaluate the dereplication steps"""
    input:
        #rules.prokaryotes__annotate__quast.input,
        rules.prokaryotes__annotate__checkm2.input,
        rules.prokaryotes__annotate__dram.input,
