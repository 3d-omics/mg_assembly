include: "functions.smk"
include: "concoct.smk"
include: "metabat2.smk"
include: "maxbin2.smk"
include: "magscot.smk"
include: "quast.smk"


rule bin_run:
    """Run the binning steps skipping evaluation"""
    input:
        rules.bin_magscot.input,


rule bin_eval:
    """Run the binning evaluation steps"""
    input:
        rules.bin_quast.input,


rule bin:
    """Run all the binning evaluation steps"""
    input:
        rules.bin_magscot.input,
        rules.bin_quast.input,
