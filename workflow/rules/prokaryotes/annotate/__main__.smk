include: "quast.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "checkm2.smk"


rule annotate:
    """Evaluate the dereplication steps"""
    input:
        rules.prokaryotes__annotate__quast.output,
        rules.prokaryotes__annotate__checkm2.output,
        rules.prokaryotes__annotate__dram.input,
