include: "functions.smk"
# include: "vamb/run.smk"
include: "concoct/run.smk"
include: "metabat2/run.smk"
include: "maxbin2/run.smk"


# include: "metawrap/run.smk"


rule bin_run:
    input:
        rules.concoct.input,
        rules.metabat2.input,
        rules.maxbin2.input,


# rule bin_eval:
#     input:
#         None


rule bin:
    input:
        rules.bin_run.input,
