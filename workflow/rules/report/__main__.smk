include: "__functions__.smk"
include: "step.smk"

rule report:
    """Report by step and by assembly"""
    input:
        rules.report__step.input,
