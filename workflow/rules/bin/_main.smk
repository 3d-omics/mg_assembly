include: "_functions.smk"
include: "concoct.smk"
include: "metabat2.smk"
include: "maxbin2.smk"
# inlcude: "metawrap2.smk"
include: "magscot.smk"
include: "quast.smk"


rule bin__run:
    """Run the binning steps skipping evaluation"""
    input:
        rules.bin__magscot.input,


rule bin__eval:
    """Run the binning evaluation steps"""
    input:
        rules.bin__quast.input,


rule bin:
    """Run all the binning evaluation steps"""
    input:
        rules.bin__magscot.input,
        rules.bin__quast.input,
