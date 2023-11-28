rule _dereplicate__drep__separate_bins:
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


rule _dereplicate__drep__run:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=DREP / "separated_bins",
    output:
        dereplicated_genomes=directory(DREP / "dereplicated_genomes"),
        data=DREP / "data.tar.gz",
        data_tables=DREP / "data_tables.tar.gz",
    log:
        DREP / "drep.log",
    conda:
        "drep.yml"
    threads: 24
    params:
        out_dir=DREP,
    resources:
        mem_mb=64 * 1024,
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

        for folder in data data_tables ; do
            tar \
                --create \
                --directory {params.out_dir} \
                --file {params.out_dir}/${{folder}}.tar.gz \
                --remove-files \
                --use-compress-program="pigz --processes {threads}" \
                --verbose \
                ${{folder}} \
            2>> {log} 1>&2
        done
        """


rule _dereplicate__drep__join_genomes:
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


rule dereplicate__drep:
    input:
        DREP / "dereplicated_genomes.fa",
