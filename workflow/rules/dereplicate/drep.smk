rule dereplicate_drep:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=[MAGSCOT / f"{assembly_id}.fa" for assembly_id in ASSEMBLIES],
    output:
        out_dir=directory(DREP / "dereplicated_genomes"),
    log:
        DEREPLICATE / "drep.log",
    conda:
        "dereplicate.yml"
    threads: params["dereplicate"]["drep"]["threads"]
    params:
        out_dir=DREP,
    resources:
        mem_mb=double_ram(params["dereplicate"]["drep"]["memory_gb"]),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        dRep dereplicate \
            {params.out_dir} \
            --processors {threads} \
            --completeness 50 \
            --S_ani 0.9 \
            --genomes {input.genomes} \
        2> {log}
        """


rule dereplicate_join_genomes:
    """Join all the dereplicated genomes into a single file."""
    input:
        DREP / "dereplicated_genomes",
    output:
        DREP / "dereplicated_genomes.fa",
    log:
        DREP / "dereplicated_genomes.log",
    conda:
        "dereplicate.yml"
    threads: 1
    shell:
        """
        cat {input}/*.fa > {output} 2> {log}
        """
