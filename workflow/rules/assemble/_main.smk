include: "_functions.smk"
include: "megahit.smk"
include: "bowtie2.smk"
include: "samtools.smk"
include: "coverm.smk"
include: "quast.smk"


rule assemble__run:
    """Run all assembly rules (no evaluation)"""
    input:
        rules.assemble__bowtie2.input,


rule assemble__eval:
    """Run the assemble evaluation rules"""
    input:
        rules.assemble__quast.input,
        rules.assemble__coverm.input,
        rules.assemble__samtools.input,


rule assemble:
    """Run the assemble module"""
    input:
        rules.assemble__eval.input,
        rules.assemble__run.input,
