rule dereplicate_drep_separate_bins:
    input:
        assemblies=[MAGSCOT / f"{assembly_id}.fa" for assembly_id in ASSEMBLIES],
    output:
        out_dir=directory(DREP / "separated_bins"),
    log:
        DREP / "separate_bins.log",
    conda:
        "drep.yml"
    shell:
        """
        mkdir --parents {output.out_dir} 2> {log} 1>&2

        ( cat {input.assemblies} \
        | paste - - \
        | tr -d ">" \
        | tr "@" "\t" \
        | awk \
            '{{print ">" $1 "@" $2 "\\n" $3 > "{output.out_dir}/" $1 ".fa" }}' \
        ) >> {log} 2>&1
        """


rule dereplicate_drep_run:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=DREP / "separated_bins",
    output:
        dereplicated_genomes=directory(DREP / "dereplicated_genomes"),
    log:
        DREP / "drep.log",
    conda:
        "drep.yml"
    threads: params["dereplicate"]["drep"]["threads"]
    params:
        out_dir=DREP,
    resources:
        mem_mb=double_ram(params["dereplicate"]["drep"]["memory_gb"]),
        runtime=7 * 24 * 60,
    retries: 5
    shell:
        """
        dRep dereplicate \
            {params.out_dir} \
            --processors {threads} \
            --completeness 50 \
            --S_ani 0.9 \
            --genomes {input.genomes}/*.fa \
        2>> {log} 1>&2

        tar \
            --create \
            --verbose \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --file \
            {params.out_dir}/data.tar.gz \
            {params.out_dir}/data \
        2>> {log} 1>&2

        tar \
            --create \
            --verbose \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --file \
            {params.out_dir}/data_tables.tar.gz \
            {params.out_dir}/data_tables \
        2>> {log}
        """


rule dereplicate_drep_join_genomes:
    """Join all the dereplicated genomes into a single file."""
    input:
        DREP / "dereplicated_genomes",
    output:
        DREP / "dereplicated_genomes.fa",
    log:
        DREP / "dereplicated_genomes.log",
    conda:
        "drep.yml"
    threads: 1
    shell:
        """
        cat {input}/*.fa > {output} 2> {log}
        """


rule dereplicate_drep:
    input:
        DREP / "dereplicated_genomes.fa",
