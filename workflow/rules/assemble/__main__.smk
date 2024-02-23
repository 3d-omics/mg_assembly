include: "__functions__.smk"
include: "megahit.smk"


rule assemble:
    input:
        rules.assemble__megahit.input,
