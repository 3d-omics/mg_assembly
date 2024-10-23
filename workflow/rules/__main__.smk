include: "folders.smk"
include: "__functions__.smk"
include: "assemble/__main__.smk"
include: "prokaryotes/__main__.smk"
include: "viruses/__main__.smk"


module preprocess:
    snakefile:
        github("jlanga/mg_preprocess", path="workflow/Snakefile", branch="devel")
    config:
        params


use rule * from preprocess
