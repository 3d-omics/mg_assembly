include: "functions.smk"
include: "concoct.smk"
include: "metabat2.smk"
include: "maxbin2.smk"
include: "magscot.smk"
include: "quast.smk"


rule bin_run:
    input:
        rules.bin_magscot.input,


rule bin_eval:
    input:
        rules.bin_quast.input,


rule bin:
    input:
        rules.bin_magscot.input,
        rules.bin_quast.input,
