include: "__functions__.smk"
include: "megahit.smk"


rule assemble:
    """Run everything in the assemble module"""
    input:
        rules.assemble__megahit.input,
