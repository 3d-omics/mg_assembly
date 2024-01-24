rule _viral__htseq__count_genes:
    input:
        crams=get_crams_from_assembly_id,
        genes_gff=GENOMAD / "{assembly_id}_summary" / "{assembly_id}_virus_genes.gff",
    output:
        tsv=HTSEQ / "{assembly_id}_genes.tsv",
    log:
        HTSEQ / "{assembly_id}_genes.log",
    conda:
        "__environment__.yml"
    shell:
        """
        htseq-count \
            --order pos \
            --type CDS \
            --idattr ID \
            --counts_output {output.tsv} \
            {input.crams} \
            {input.genes_gff} \
        2> {log} 1>&2
        """


rule viral__htseq__genes:
    input:
        [HTSEQ / f"{assembly_id}_genes.tsv" for assembly_id in ASSEMBLIES],


rule _viral__htseq__contigs:
    input:
        crams=get_crams_from_assembly_id,
        contigs_gff=GENOMAD
        / "{assembly_id}_summary"
        / "{assembly_id}_virus_contigs.gff",
    output:
        tsv=HTSEQ / "{assembly_id}_contigs.tsv",
    log:
        HTSEQ / "{assembly_id}_contigs.log",
    conda:
        "__environment__.yml"
    shell:
        """
        htseq-count \
            --order pos \
            --type CDS \
            --idattr ID \
            --counts_output {output.tsv} \
            {input.crams} \
            {input.contigs_gff} \
        2> {log} 1>&2
        """


rule viral__htseq__contigs:
    input:
        [HTSEQ / f"{assembly_id}_contigs.tsv" for assembly_id in ASSEMBLIES],


rule viral__htseq:
    input:
        rules.viral__htseq__genes.input,
        rules.viral__htseq__contigs.input,
