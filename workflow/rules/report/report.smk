include: "functions.smk"
include: "by_step.smk"
include: "by_assembly.smk"


rule report:
    input:
        rules.report_step.input,
        rules.report_assembly.input,
