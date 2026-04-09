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


use rule concatenate__flat_to_gzipped as viruses__cluster__genomad__concatenate_plasmid_fna with:
    input:
        [VIR_GENOMADC / f"{assembly_id}_plasmid.fna" for assembly_id in ASSEMBLIES]
        + ["/dev/null"],
    output:
        VIR_CLUSTER / "genomad_plasmid.fna.gz",
    log:
        VIR_CLUSTER / "genomad_plasmid.log",


use rule concatenate__flat_to_gzipped as viruses__cluster__genomad__concatenate_plasmid_proteins_faa with:
    input:
        [
            VIR_GENOMADC / f"{assembly_id}_plasmid_proteins.faa"
            for assembly_id in ASSEMBLIES
        ]
        + ["/dev/null"],
    output:
        VIR_CLUSTER / "genomad_plasmid_proteins.faa.gz",
    log:
        VIR_CLUSTER / "genomad_plasmid_proteins.log",


use rule concatenate__flat_to_gzipped as viruses__cluster__genomad__concatenate_virus_fna with:
    input:
        [VIR_GENOMADC / f"{assembly_id}_virus.fna" for assembly_id in ASSEMBLIES]
        + ["/dev/null"],
    output:
        VIR_CLUSTER / "genomad_virus.fna.gz",
    log:
        VIR_CLUSTER / "genomad_virus.log",


use rule concatenate__flat_to_gzipped as viruses__cluster__genomad__concatenate_virus_proteins_faa with:
    input:
        [
            VIR_GENOMADC / f"{assembly_id}_virus_proteins.faa"
            for assembly_id in ASSEMBLIES
        ]
        + ["/dev/null"],
    output:
        VIR_CLUSTER / "genomad_virus_proteins.faa.gz",
    log:
        VIR_CLUSTER / "genomad_virus_proteins.log",


rule viruses__cluster__genomad__concatenate_fastas:
    input:
        VIR_CLUSTER / "genomad_plasmid.fna.gz",
        VIR_CLUSTER / "genomad_plasmid_proteins.faa.gz",
        VIR_CLUSTER / "genomad_virus.fna.gz",
        VIR_CLUSTER / "genomad_virus_proteins.faa.gz",


use rule csvtk__concat as viruses__cluster__genomad__concatenate__plasmid_genes with:
    input:
        [
            VIR_GENOMADC / f"{assembly_id}_plasmid_genes.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + ["/dev/null"],
    output:
        VIR_CLUSTER / "genomad_plasmid_genes.tsv.gz",
    log:
        VIR_CLUSTER / "genomad_plasmid_genes.log",


use rule csvtk__concat as viruses__cluster__genomad__concatenate__plasmid_summary with:
    input:
        [
            VIR_GENOMADC / f"{assembly_id}_plasmid_summary.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + ["/dev/null"],
    output:
        VIR_CLUSTER / "genomad_plasmid_summary.tsv.gz",
    log:
        VIR_CLUSTER / "genomad_plasmid_summary.log",


use rule csvtk__concat as viruses__cluster__genomad__concatenate__virus_genes with:
    input:
        [VIR_GENOMADC / f"{assembly_id}_virus_genes.tsv" for assembly_id in ASSEMBLIES]
        + ["/dev/null"],
    output:
        VIR_CLUSTER / "genomad_virus_genes.tsv.gz",
    log:
        VIR_CLUSTER / "genomad_virus_genes.log",


use rule csvtk__concat as viruses__cluster__genomad__concatenate__virus_summary with:
    input:
        [
            VIR_GENOMADC / f"{assembly_id}_virus_summary.tsv"
            for assembly_id in ASSEMBLIES
        ]
        + ["/dev/null"],
    output:
        VIR_CLUSTER / "genomad_virus_summary.tsv.gz",
    log:
        VIR_CLUSTER / "genomad_virus_summary.log",


rule viruses__cluster__genomad__concatenate_tsvs:
    input:
        VIR_CLUSTER / "genomad_plasmid_genes.tsv.gz",
        VIR_CLUSTER / "genomad_plasmid_summary.tsv.gz",
        VIR_CLUSTER / "genomad_virus_genes.tsv.gz",
        VIR_CLUSTER / "genomad_virus_summary.tsv.gz",


rule viruses__cluster__genomad__all:
    input:
        rules.viruses__cluster__genomad__concatenate_fastas.output,
        rules.viruses__cluster__genomad__concatenate_tsvs.output,
