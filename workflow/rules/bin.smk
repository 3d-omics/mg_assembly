include: "bin/vamb.smk"
include: "bin/concoct.smk"
include: "bin/metabat2.smk"
include: "bin/maxbin2.smk"
include: "bin/metawrap.smk"


rule bin:
    input:
        rules.bin_metawrap.output,
