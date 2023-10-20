include: "coverm.smk"
include: "quast.smk"
include: "samtools.smk"


rule assemble_eval:
    input:
        rules.assemble_eval_quast.input,
        rules.assemble_eval_coverm_contig.input,
        rules.assemble_eval_coverm_genome.input,
        rules.assemble_eval_samtools.input,
