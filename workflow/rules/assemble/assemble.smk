include: "functions.smk"
include: "megahit.smk"
include: "bowtie2.smk"
include: "samtools.smk"
include: "coverm.smk"
include: "quast.smk"


rule assemble_run:
    """Run all assembly rules (no evaluation)"""
    input:
        rules.assemble_bowtie2.input,


rule assemble_eval:
    """Run the assemble evaluation rules"""
    input:
        rules.assemble_quast.input,
        rules.assemble_coverm.input,
        rules.assemble_samtools.input,


rule assemble:
    """Run the assemble module"""
    input:
        rules.assemble_eval.input,
        rules.assemble_run.input,
