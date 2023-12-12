include: "__functions__.smk"
include: "quast.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "checkm2.smk"


rule annotate:
    """Evaluate the dereplication steps"""
    input:
        rules.annotate__quast.output,
        rules.annotate__checkm2.output,


rule annotate__with_gtdbtk:
    """Run the evaluation steps + GTDB-Tk"""
    input:
        rules.annotate.input,
        rules.annotate__gtdbtk.output,


rule annotate__with_dram:
    """Run the evaluation steps + DRAM"""
    input:
        rules.annotate.input,
        rules.annotate__dram.input,
