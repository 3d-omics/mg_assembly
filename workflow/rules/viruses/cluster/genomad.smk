rule viruses__cluster__genomad__download_database:
    output:
        directory(features["databases"]["genomad"]),
    log:
        f"{features["databases"]["genomad"]}.log",
    conda:
        ENVS / "genomad.yml"
    shell:
        """
        genomad download-database \
            --verbose \
            $(dirname {output}) \
        2> {log} 1>&2

        mv \
            --verbose \
            $(dirname {output})/genomad_db \
            {output} \
        2>> {log} 1>&2
        """


rule viruses__cluster__genomad__run:
    input:
        fasta=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
        database=features["databases"]["genomad"],
    output:
        plasmid=temp(VIR_GENOMADC / "{assembly_id}_plasmid.fna"),
        plasmid_genes=temp(VIR_GENOMADC / "{assembly_id}_plasmid_genes.tsv"),
        plasmid_proteins=temp(VIR_GENOMADC / "{assembly_id}_plasmid_proteins.faa"),
        plasmid_summary=temp(VIR_GENOMADC / "{assembly_id}_plasmid_summary.tsv"),
        json=VIR_GENOMADC / "{assembly_id}_summary.json",
        virus=temp(VIR_GENOMADC / "{assembly_id}_virus.fna"),
        virus_genes=temp(VIR_GENOMADC / "{assembly_id}_virus_genes.tsv"),
        virus_proteins=temp(VIR_GENOMADC / "{assembly_id}_virus_proteins.faa"),
        virus_summary_tsv=temp(VIR_GENOMADC / "{assembly_id}_virus_summary.tsv"),
    log:
        VIR_GENOMADC / "{assembly_id}.log",
    conda:
        ENVS / "genomad.yml"
    params:
        filtering=params["viral"]["genomad"]["filtering"],
        genomad_workdir=VIR_GENOMADC,
        genomad_summary_dir=lambda w: VIR_GENOMADC / f"{w.assembly_id}_summary",
        extra=params["viral"]["genomad"]["extra"],
        use_cuda=params["viral"]["genomad"]["use_cuda"],
    shadow:
        "minimal"
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
    shell:
        """
        if [[ $(gzip -dc {input.fasta} | wc -l ) -lt 2 ]] ; then
            echo "Empty fasta. Touching outputs" 2> {log}  1>&2
            touch {output}
            exit 0
        fi

        {params.use_cuda}

        genomad end-to-end \
            {params.filtering} \
            --cleanup \
            --restart \
            --verbose \
            --threads {threads} \
            {params.extra} \
            {input.fasta} \
            {params.genomad_workdir} \
            {input.database} \
        2> {log} 1>&2

        mv \
            --verbose \
            {params.genomad_summary_dir}/* \
            {params.genomad_workdir} \
        2>> {log} 1>&2
        """


rule viruses__cluster__genomad__concatenate_fastas:
    input:
        plasmid_fnas=[
            VIR_GENOMADC / f"{assembly_id}_plasmid.fna" for assembly_id in ASSEMBLIES
        ],
        plasmid_proteins_fnas=[
            VIR_GENOMADC / f"{assembly_id}_plasmid_proteins.faa"
            for assembly_id in ASSEMBLIES
        ],
        virus_fnas=[
            VIR_GENOMADC / f"{assembly_id}_virus.fna" for assembly_id in ASSEMBLIES
        ],
        virus_proteins_fnas=[
            VIR_GENOMADC / f"{assembly_id}_virus_proteins.faa"
            for assembly_id in ASSEMBLIES
        ],
    output:
        plasmid_fna=VIR_CLUSTER / "genomad_plasmid.fna.gz",
        plasmid_proteins_fna=VIR_CLUSTER / "genomad_plasmid_proteins.fna.gz",
        virus_fna=VIR_CLUSTER / "genomad_virus.fna.gz",
        virus_proteins_fna=VIR_CLUSTER / "genomad_virus_proteins.fna.gz",
    log:
        VIR_CLUSTER / "genomad.fastas.log",
    conda:
        ENVS / "genomad.yml"
    threads: 24
    shell:
        """
        (
            cat {input.plasmid_fnas} | bgzip --threads {threads} > {output.plasmid_fna}
            cat {input.plasmid_proteins_fnas} | bgzip --threads {threads} > {output.plasmid_proteins_fna}
            cat {input.virus_fnas} | bgzip --threads {threads} > {output.virus_fna}
            cat {input.virus_proteins_fnas} | bgzip --threads {threads} > {output.virus_proteins_fna}
        ) 2> {log}
        """


rule viruses__cluster__genomad__aggregate_tsvs:
    input:
        plasmid_genes=[
            VIR_GENOMADC / f"{assembly_id}_plasmid_genes.tsv"
            for assembly_id in ASSEMBLIES
        ],
        plasmid_summary=[
            VIR_GENOMADC / f"{assembly_id}_plasmid_summary.tsv"
            for assembly_id in ASSEMBLIES
        ],
        virus_genes=[
            VIR_GENOMADC / f"{assembly_id}_virus_genes.tsv"
            for assembly_id in ASSEMBLIES
        ],
        virus_summary_tsv=[
            VIR_GENOMADC / f"{assembly_id}_virus_summary.tsv"
            for assembly_id in ASSEMBLIES
        ],
    output:
        plasmid_genes=VIR_CLUSTER / "genomad_plasmid_genes.tsv.gz",
        plasmid_summary=VIR_CLUSTER / "genomad.plasmid_summary.tsv.gz",
        virus_genes=VIR_CLUSTER / "genomad.virus_genes.tsv.gz",
        virus_summary_tsv=VIR_CLUSTER / "genomad.virus_summary.tsv.gz",
    log:
        VIR_CLUSTER / "genomad.tsvs.log",
    threads: 24
    conda:
        ENVS / "genomad.yml"
    shell:
        """
        csvtk concat \
            --tabs \
            {input.plasmid_genes} \
        | bgzip --compress-level 0 --threads {threads} \
        > {output.plasmid_genes}

        csvtk concat \
            --tabs \
            {input.plasmid_summary} \
        | bgzip --compress-level 0 --threads {threads} \
        > {output.plasmid_summary}

        csvtk concat \
            --tabs \
            {input.virus_genes} \
        | bgzip --compress-level 0 --threads {threads} \
        > {output.virus_genes}

        csvtk concat \
            --tabs \
            {input.virus_summary_tsv} \
        | bgzip --compress-level 0 --threads {threads} \
        > {output.virus_summary_tsv}
        """


rule viruses__cluster__genomad__all:
    input:
        rules.viruses__cluster__genomad__concatenate_fastas.output,
        rules.viruses__cluster__genomad__aggregate_tsvs.output,
