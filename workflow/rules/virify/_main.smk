include: "virify.yml"


rule virify:
    """Run virify"""
    input:
        rules.virify_all.input,
