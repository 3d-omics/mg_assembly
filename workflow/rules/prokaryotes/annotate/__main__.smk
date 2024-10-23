include: "mags.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "checkm2.smk"
include: "drep.smk"
include: "quast.smk"


rule prokaryotes__annotate__all:
    """Evaluate the dereplication steps"""
    input:
        rules.prokaryotes__annotate__checkm2__all.input,
        rules.prokaryotes__annotate__dram__all.input,
        rules.prokaryotes__annotate__gtdbtk__all.input,
        rules.prokaryotes__annotate__drep__all.input,
