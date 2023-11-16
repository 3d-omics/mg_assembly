include: "_functions.smk"
include: "drep.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "quast.smk"
include: "samtools.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "checkm2.smk"


rule dereplicate_run:
    """Run the dereplication steps without evaluation"""
    input:
        rules.dereplicate_bowtie2.input,


rule dereplicate_eval:
    """Evaluate the dereplication steps"""
    input:
        rules.dereplicate_coverm.input,
        rules.dereplicate_samtools.input,
        rules.dereplicate_quast.output,
        rules.dereplicate_checkm2.output,


rule dereplicate_eval_with_dram:
    """Run the evaluation steps + DRAM"""
    input:
        rules.dereplicate_eval.input,
        rules.dereplicate_dram.input,


rule dereplicate_eval_with_gtdbtk:
    """Run the evaluation steps + GTDB-Tk"""
    input:
        rules.dereplicate_eval.input,
        rules.dereplicate_gtdbtk.output,


rule dereplicate:
    """Run the dereplication steps with evaluation"""
    input:
        rules.dereplicate_run.input,
        rules.dereplicate_eval.input,
