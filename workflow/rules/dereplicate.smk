rule dereplicate_compose_scores:
    input:
        tables=[MAGSCOT / "{assembly_id}/magscot.scores.out"],
    output:
        table=DEREPLICATE / "scores.tsv",
    log:
        DEREPLICATE / "scores.log",
    conda:
        "../envs/dereplicate.yml"
    shell:
        """
        """


rule dereplicate_drep:
    input:
        genomes=[MAGSCOT / f"{assembly_id}.fa" for assembly_id in ASSEMBLIES],
    output:
        out_dir=directory(DREP),
        genome_info=DEREPLICATE / "genome_info.csv",
    log:
        DEREPLICATE / "drep.log",
    conda:
        "../envs/dereplicate.yml"
    shell:
        """
        dRep dereplicate \
            --processors {threads} \
            --completeness 50 \
            --S_ani 0.9 \
            --genomes {input.genomes} \
            --genomeInfo {output.genome_info} \
            {output.out_dir} \
        2> {log}
        """


rule dereplicate:
    input:
        DEREPLICATE / "genome_info.csv",
