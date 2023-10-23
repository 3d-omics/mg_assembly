include: "functions.smk"
include: "concoct.smk"
include: "metabat2.smk"
include: "maxbin2.smk"
include: "magscot.smk"


rule bin_run:
    input:
        rules.bin_magscot.input,


# rule bin_eval:
#     input:
#         None


rule bin:
    input:
        rules.bin_magscot.input,
