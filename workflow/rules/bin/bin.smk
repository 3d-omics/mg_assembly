include: "functions.smk"
include: "concoct/concoct.smk"
include: "metabat2/metabat2.smk"
include: "maxbin2/maxbin2.smk"
include: "magscot/magscot.smk"


rule bin_run:
    input:
        rules.magscot.input,


# rule bin_eval:
#     input:
#         None


rule bin:
    input:
        rules.magscot.input,
