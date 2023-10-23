include: "functions.smk"
include: "by_step.smk"
include: "by_assembly.smk"


rule report:
    """Report by step and by assembly"""
    input:
        rules.report_step.input,
        rules.report_assembly.input,
