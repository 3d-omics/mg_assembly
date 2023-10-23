include: "coverm.smk"
include: "quast.smk"
include: "samtools.smk"
include: "gtdbtk.smk"
include: "dram.smk"


rule dereplicate_eval:
    input:
        rules.dereplicate_coverm_genome.output,
        rules.dereplicate_coverm_contig.output,
        rules.dereplicate_samtools.input,
        rules.dereplicate_quast.output,


rule dereplicate_eval_with_dram:
    input:
        rules.dereplicate_eval.input,
        rules.dereplicate_dram_distill.output,


rule dereplicate_eval_with_gtdbtk:
    input:
        rules.dereplicate_eval.input,
        rules.dereplicate_gtdbtk.output,
