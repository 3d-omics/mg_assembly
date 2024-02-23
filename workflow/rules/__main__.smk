include: "folders.smk"
include: "__functions__.smk"
include: "reads/__main__.smk"
include: "reference/__main__.smk"
include: "preprocess/__main__.smk"
include: "assemble/__main__.smk"
include: "prokaryotes/__main__.smk"
include: "viruses/__main__.smk"
include: "report/__main__.smk"

module helpers:
    snakefile: "helpers/Snakefile"

use rule * from helpers as helpers__*
