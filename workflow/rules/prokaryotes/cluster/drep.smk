rule prokaryotes__cluster__drep__separate_bins__:
    input:
        assemblies=[MAGSCOT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
    output:
        out_dir=directory(DREP / "separated_bins"),
    log:
        DREP / "separate_bins.log",
    conda:
        "__environment__.yml"
    threads: 24
    shell:
        """
        mkdir --parents {output.out_dir} 2> {log} 1>&2

        ( gzip \
            --decompress \
            --stdout \
            {input.assemblies} \
        | paste - - \
        | tr -d ">" \
        | tr "@" "\t" \
        | awk \
            '{{print ">" $1 "@" $2 "\\n" $3 > "{output.out_dir}/" $1 ".fa" }}' \
        ) >> {log} 2>&1

        pigz \
            --best \
            --processes {threads} \
            --verbose \
            {output.out_dir}/*.fa \
        2>> {log} 1>&2
        """


rule prokaryotes__cluster__drep__run__:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=DREP / "separated_bins",
    output:
        fasta=DREP / "dereplicated_genomes.fa.gz",
        dereplicated_genomes=directory(DREP / "dereplicated_genomes"),
        data=DREP / "data.tar.gz",
        data_tables=DREP / "data_tables.tar.gz",
    log:
        DREP / "drep.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        out_dir=DREP,
    resources:
        mem_mb=64 * 1024,
        runtime=7 * 24 * 60,
        attempt=get_attempt,
    retries: 5
    shell:
        """
        rm \
            --recursive \
            --force \
            {params.out_dir}/data_tables \
            {params.out_dir}/data \
            {params.out_dir}/dereplicated_genomes \
            {params.out_dir}/figures \
            {params.out_dir}/log \
        2> {log}.{resources.attempt} 1>&2

        pigz \
            --decompress \
            --keep \
            {input.genomes}/*.fa.gz \
        2>> {log}.{resources.attempt} 1>&2

        dRep dereplicate \
            {params.out_dir} \
            --processors {threads} \
            --completeness 50 \
            --S_ani 0.9 \
            --genomes {input.genomes}/*.fa \
        2>> {log}.{resources.attempt} 1>&2

        ( cat \
            {params.out_dir}/dereplicated_genomes/*.fa \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.fasta} \
        ) 2>> {log}.{resources.attempt} 1>&2

        bgzip \
            --level 9 \
            --processes {threads} \
            --verbose \
            {params.out_dir}/dereplicated_genomes/*.fa \
        2>> {log}.{resources.attempt} 1>&2

        rm \
            --force \
            --verbose \
            {input.genomes}/*.fa \
        2>> {log}.{resources.attempt} 1>&2

        for folder in data data_tables ; do
            tar \
                --create \
                --directory {params.out_dir} \
                --file {params.out_dir}/${{folder}}.tar.gz \
                --remove-files \
                --use-compress-program="pigz --processes {threads}" \
                --verbose \
                ${{folder}} \
            2>> {log}.{resources.attempt} 1>&2
        done

        mv {log}.{resources.attempt} {log}
        """


rule prokaryotes__cluster__drep:
    input:
        DREP / "dereplicated_genomes.fa.gz",
