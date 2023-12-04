include: "_functions.smk"
include: "quast.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "checkm2.smk"


rule dereplicate__eval:
    """Evaluate the dereplication steps"""
    input:
        rules.dereplicate__quast.output,
        rules.dereplicate__checkm2.output,


rule dereplicate__eval_with_gtdbtk:
    """Run the evaluation steps + GTDB-Tk"""
    input:
        rules.dereplicate__eval.input,
        rules.dereplicate__gtdbtk.output,


rule dereplicate__eval_with_dram:
    """Run the evaluation steps + DRAM"""
    input:
        rules.dereplicate__eval.input,
        rules.dereplicate__dram.input,


rule dereplicate:
    """Run the dereplication steps with evaluation"""
    input:
        rules.dereplicate__eval.input,
