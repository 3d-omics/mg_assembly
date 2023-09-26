include: "report/functions.smk"
include: "report/by_step.smk"
include: "report/by_assembly.smk"


rule report:
    input:
        rules.report_step.input,
        rules.report_assembly.input,
