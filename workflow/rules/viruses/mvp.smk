include: "mvp/functions.smk"
include: "mvp/input_prep.smk"
include: "mvp/00_set_up.smk"
include: "mvp/01_genomad_checkv.smk"
include: "mvp/02_filter.smk"
include: "mvp/03_clustering.smk"
include: "mvp/04_read_mapping.smk"
include: "mvp/05_votu_table.smk"
include: "mvp/06_functional_annotation.smk"
include: "mvp/07_binning.smk"


rule viruses__mvp__all:
    input:
        rules.viruses__mvp__01_genomad_checkv__all.input,
        rules.viruses__mvp__02_filter__all.input,
        rules.viruses__mvp__03_clustering.output,
        rules.viruses__mvp__04_read_mapping__all.input,
        rules.viruses__mvp__05_votu_table.output,
        rules.viruses__mvp__06_functional_annotation.output,
        rules.viruses__mvp__07_binning.output,
