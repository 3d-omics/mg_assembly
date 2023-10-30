rule dereplicate_drep_separate_bins:
    input:
        assemblies=[MAGSCOT / f"{assembly_id}.fa" for assembly_id in ASSEMBLIES],
    output:
        out_dir=directory(DREP / "separated_bins"),
    log:
        DEREPLICATE / "separate_bins.log",
    conda:
        "dereplicate.yml"
    shell:
        """
        cat /dev/null > {log}

        ( cat {input.assemblies} \
        | paste - - \
        | tr -d ">" \
        | tr "@:" "\t" \
        | awk \
            -v out_dir="{output.out_dir}" \
            '{{print ">" $1 ":" $2 "@" $3 "\n" $4 > out_dir "/" $1 ":" $2 ".fa" }}' \
        ) >> {log} 2>&1
        """


rule dereplicate_drep:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=DREP / "separated_bins",
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
            {output.out_dir} \
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
