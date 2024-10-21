include: "folders.smk"
include: "__functions__.smk"
include: "assemble/__main__.smk"
include: "prokaryotes/__main__.smk"
include: "viruses/__main__.smk"
include: "report/__main__.smk"


module helpers:
    snakefile:
        "helpers/Snakefile"

module preprocess:
    snakefile:
        github("jlanga/mg_preprocess", path="workflow/Snakefile", branch="devel")
    config:
        params


use rule * from preprocess


use rule * from helpers as helpers__*
