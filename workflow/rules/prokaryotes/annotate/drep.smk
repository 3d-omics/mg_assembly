rule prokaryotes__annotate__drep__quality_report__:
    input:
        PROK_ANN / "checkm2.quality_report.tsv",
    output:
        temp(PROK_ANN / "drep.quality_report.tsv"),
    log:
        PROK_ANN / "drep.quality_report.log",
    conda:
        "__environment__.yml"
    shell:
        """
        echo \
            -E "genome,completeness,contamination" \
        > {output} \
        2> {log}

        ( tail -n+2 {input} \
        | cut -f 1-3 \
        | awk \
            '{{print $1 ".fa," $2 "," $3}}' \
        ) \
        >> {output} \
        2>> {log}
        """


rule prokaryotes__annotate__drep__dereplicate__:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=MAGS,
        quality_report=PROK_ANN / "drep.quality_report.tsv",
    output:
        work_dir=temp(directory(PROK_ANN / "drep.{secondary_ani}.dir")),
    log:
        PROK_ANN / "drep.{secondary_ani}.log",
    conda:
        "__environment__.yml"
    params:
        secondary_ani=lambda w: w.secondary_ani,
        minimum_completeness=params["prokaryotes"]["annotate"]["drep"][
            "minimum_completeness"
        ],
        maximum_contamination=params["prokaryotes"]["annotate"]["drep"][
            "maximum_contamination"
        ],
    shell:
        """
        dRep dereplicate \
            {output.work_dir} \
            --S_ani         {params.secondary_ani} \
            --completeness  {params.minimum_completeness} \
            --contamination {params.maximum_contamination} \
            --genomes       {input.genomes}/*.fa \
            --processors    {threads} \
            --genomeInfo    {input.quality_report} \
        2>> {log} 1>&2
        """


rule prokaryotes__annotate__drep__get_fasta__:
    input:
        work_dir=PROK_ANN / "drep.{secondary_ani}.dir",
    output:
        fasta=PROK_ANN / "drep.{secondary_ani}.fa.gz",
    log:
        PROK_ANN / "drep.{secondary_ani}.fa.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( cat \
            {input.work_dir}/dereplicated_genomes/*.fa \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.fasta} \
        ) 2> {log}
        """


rule prokaryotes__annotate__drep__tarball__:
    input:
        work_dir=PROK_ANN / "drep.{secondary_ani}.dir",
    output:
        tarball=PROK_ANN / "drep.{secondary_ani}.tar.gz",
    log:
        PROK_ANN / "drep.{secondary_ani}.tar.log",
    conda:
        "__environment__.yml"
    shell:
        """
        tar \
            --create \
            --file {output.tarball} \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            {input.work_dir} \
        2>> {log} 1>&2
        """


rule prokaryotes__annotate__drep:
    input:
        [PROK_ANN / f"drep.{secondary_ani}.tar.gz" for secondary_ani in SECONDARY_ANIS],
        [PROK_ANN / f"drep.{secondary_ani}.fa.gz" for secondary_ani in SECONDARY_ANIS],


localrules:
    prokaryotes__annotate__drep__quality_report__,
