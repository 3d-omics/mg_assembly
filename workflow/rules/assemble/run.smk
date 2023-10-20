include: "megahit.smk"
include: "bowtie2.smk"


rule assemble_run:
    """Run all assembly rules (no evaluation)"""
    input:
        rules.assemble_bowtie2_all.input,
