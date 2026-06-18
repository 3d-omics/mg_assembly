include: "annotate/mags.smk"
include: "annotate/gtdbtk.smk"
include: "annotate/dram.smk"
include: "annotate/checkm2.smk"
include: "annotate/drep.smk"
include: "annotate/quast.smk"
include: "annotate/bakta.smk"


rule prokaryotes__annotate__all:
    """Evaluate the dereplication steps"""
    input:
        rules.prokaryotes__annotate__mags__all.input,
        rules.prokaryotes__annotate__checkm2__all.input,
        rules.prokaryotes__annotate__dram__all.input,
        rules.prokaryotes__annotate__gtdbtk__all.input,
        rules.prokaryotes__annotate__drep__all.input,
        rules.prokaryotes__annotate__bakta__all.input,
        rules.prokaryotes__annotate__quast__all.input,
