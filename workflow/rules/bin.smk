include: "bin/functions.smk"
# include: "bin/vamb/run.smk"
# include: "bin/concoct/run.smk"
# include: "bin/metabat2/run.smk"
# include: "bin/maxbin2/run.smk"
include: "bin/metawrap/run.smk"


rule bin_run:
    input:
        rules.bin_metawrap.input,


# rule bin_eval:
#     input:
#         None


rule bin:
    input:
        rules.bin_run.input,
