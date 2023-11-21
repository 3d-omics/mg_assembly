include: "_functions.smk"
include: "step.smk"
include: "assembly.smk"


rule report:
    """Report by step and by assembly"""
    input:
        rules.report__step.input,
        rules.report__assembly.input,
