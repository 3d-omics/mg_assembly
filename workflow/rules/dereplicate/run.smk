rule dereplicate_drep:
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
        mem_mb=params["dereplicate"]["drep"]["mem_mb"],
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


rule dereplicate_bowtie2_build_one:
    """
    Index dereplicader
    """
    input:
        contigs=DREP / "dereplicated_genomes.fa",
    output:
        mock=touch(DREP_INDEX / "dereplicated_genomes"),
    log:
        DREP_INDEX / "dereplicated_genomes.log",
    conda:
        "dereplicate.yml"
    threads: 24
    params:
        extra=params["dereplicate"]["bowtie2-build"]["extra"],
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule dereplicate_bowtie2_one:
    input:
        mock=DREP_INDEX / "dereplicated_genomes",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        reference=DREP / "dereplicated_genomes.fa",
    output:
        cram=DREP_BOWTIE2 / "{sample_id}.{library_id}.cram",
    log:
        DREP_BOWTIE2 / "{sample_id}.{library_id}.log",
    conda:
        "dereplicate.yml"
    threads: 24
    params:
        extra=params["dereplicate"]["bowtie2"]["extra"],
        samtools_mem=params["dereplicate"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    resources:
        mem_mb=32 * 1024,
    shell:
        """
        (bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule dereplicate_bowtie2:
    input:
        [
            DREP_BOWTIE2 / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule dereplicate_run:
    input:
        rules.dereplicate_bowtie2.input,
