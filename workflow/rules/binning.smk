rule binning_metawrap_binning_prepare_one:
    input:
        bam=BOWTIE2_ASSEMBLY / "{assembly_id}.bam",
    output:
        forward_=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}_1.fastq",
        reverse_=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}_2.fastq",
        bwt_index=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}.fa.bwt",
        bam=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}.bam",
    log:
        METAWRAP_BINNING / "{assembly_id}.prepare.log",
    conda:
        "../envs/binning.yml"
    shell:
        """
        echo "@" > {output.forward_}
        echo "@" > {output.reverse_}
        touch {output.bwt_index}
        ln {input.bam} {output.bam}
        """


rule binning_metawrap_binning_prepare_all:
    input:
        [METAWRAP_BINNING / f"{assembly_id}" for assembly_id in samples.assembly_id],


rule binning_metawrap_binning_one:
    """Run metawrap over one assembly group
    Note: metawrap works with fastq files, but we can trick it into working by
    creating mock fastq and reference files. h/t: Raphael Eisenhofer
    Note2: metawrap is rotten. It is written in python2 and has a lot unmetabel dependencies in conda.
    Using a singularity container instead.
    """
    input:
        bam=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}.bam",
        forward_=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}_1.fastq",
        reverse_=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}_2.fastq",
        assembly=MEGAHIT / "{assembly_id}/final.contigs.fa",
    output:
        metabat2_bins=directory(METAWRAP_BINNING / "{assembly_id}/metabat2_bins"),
        maxbin2_bins=directory(METAWRAP_BINNING / "{assembly_id}/maxbin2_bins"),
        concoct_bins=directory(METAWRAP_BINNING / "{assembly_id}/concoct_bins"),
    log:
        METAWRAP_BINNING / "{assembly_id}.log",
    singularity:
        "https://depot.galaxyproject.org/singularity/metawrap-mg:1.3.0--hdfd78af_1"
    threads: 24
    params:
        min_length=params["binning"]["metawrap_binning"]["min_length"],
        extra=params["binning"]["metawrap_binning"]["extra"],
        out_folder=lambda wildcards: METAWRAP_BINNING / f"{wildcards.assembly_id}",
    shell:
        """
        metawrap binning \
            -o {params.out_folder} \
            -a {input.assembly} \
            -l {params.min_length} \
            --metabat2 \
            --maxbin2 \
            --concoct \
            {params.extra} \
            {input.forward_} \
            {input.reverse_} \
        2> {log} 1>&2
        """


rule binning_metawrap_binning_all:
    input:
        [
            METAWRAP_BINNING / f"{assembly_id}/{binner}"
            for assembly_id in samples.assembly_id
            for binner in ["concoct_bins", "maxbin2_bins", "metabat2_bins"]
        ],


rule binning_metawrap_bin_refinement_one:
    input:
        metabat2_bins=METAWRAP_BINNING / "{assembly_id}/metabat2_bins",
        maxbin2_bins=METAWRAP_BINNING / "{assembly_id}/maxbin2_bins",
        concoct_bins=METAWRAP_BINNING / "{assembly_id}/concoct_bins",
    output:
        stats=METAWRAP_REFINEMENT / "{assembly_id}_bins.stats",
        contigs=METAWRAP_REFINEMENT / "{assembly_id}_bins.contigs",
        working_folder=directory(METAWRAP_REFINEMENT / "{assembly_id}"),
        bins_folder=directory(METAWRAP_REFINEMENT / "{assembly_id}_bins"),
    log:
        METAWRAP_REFINEMENT / "{assembly_id}.log",
    singularity:
        "https://depot.galaxyproject.org/singularity/metawrap-mg:1.3.0--hdfd78af_1"
    threads: 24
    params:
        completeness=params["binning"]["metawrap_bin_refinement"]["completeness"],
        contamination=params["binning"]["metawrap_bin_refinement"]["contamination"],
        extra=params["binning"]["metawrap_bin_refinement"]["extra"],
        output_prefix=compose_metawrap_working_folder,
    resources:
        mem_mb=16 * 1024,
    shell:
        """
        metawrap bin_refinement \
            -m $(({resources.mem_mb} / 1024)) \
            -o {output.working_folder} \
            -t {threads} \
            -A {input.metabat2_bins} \
            -B {input.maxbin2_bins} \
            -C {input.concoct_bins} \
            -c {params.completeness} \
            -x {params.contamination} \
            {params.extra} \
        2> {log} 1>&2

        cp \
            {params.output_prefix}.stats \
            {output.stats} \
        2>> {log} 1>&2

        cp \
            {params.output_prefix}.contigs \
            {output.contigs} \
        2>> {log} 1>&2

        cp --recursive \
            {params.output_prefix} \
            {output.bins_folder} \
        2>> {log} 1>&2
        """


rule binning_metawrap_bin_refinement_all:
    input:
        [
            METAWRAP_REFINEMENT / f"{assembly_id}.contigs"
            for assembly_id in samples.assembly_id
        ],


rule binning_metawrap_renaming_one:
    """

    Note: doing this separatedly from the binning step because we need seqtk and it is outside the metawrap singularity container
    """
    input:
        bin_folder=METAWRAP_REFINEMENT / "{assembly_id}_bins",
    output:
        fa=METAWRAP_RENAMING / "{assembly_id}.fa",
    log:
        METAWRAP_RENAMING / "{assembly_id}.log",
    conda:
        "../envs/binning.yml"
    params:
        assembly_id=lambda wildcards: f"{wildcards.assembly_id}",
    shell:
        """
        (for bin in {input.bin_folder}/*.fa ; do
            bin_name=$(basename $bin .fa)
            setqk rename \
                {params.assembly_id}-${{bin_name}}- \
                $bin
        done > {output.fa}) 2> {log}
        """


rule binning_metawrap_renaming_all:
    input:
        [METAWRAP_RENAMING / f"{assembly_id}.fa" for assembly_id in samples.assembly_id],


rule binning_bowtie2_index_bin_one:
    input:
        METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        touch(BOWTIE2_INDEXES_BINNING / "{assembly_id}"),
    log:
        BOWTIE2_INDEXES_BINNING / "{assembly_id}.log",
    conda:
        "../envs/binning.yml"
    threads: 24
    params:
        extra=params["assembly"]["bowtie2-build"]["extra"],
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule binning_bowtie2_map_one_library_to_one_bin:
    input:
        mock=BOWTIE2_INDEXES_BINNING / "{assembly_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        reference=METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        cram=BOWTIE2_BINNING / "{assembly_id}/{sample_id}.{library_id}.cram",
    log:
        BOWTIE2_BINNING / "{assembly_id}/{sample_id}.{library_id}.log",
    conda:
        "../envs/binning.yml"
    threads: 24
    params:
        extra=params["binning"]["bowtie2"]["extra"],
        samtools_mem=params["binning"]["samtools"]["mem"],
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


rule binning_bowtie2_all:
    input:
        [
            BOWTIE2_BINNING / f"{assembly_id}/{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in samples[
                ["assembly_id", "sample_id", "library_id"]
            ].values.tolist()
        ],


rule binning_coverm_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=BOWTIE2_ASSEMBLY / "{assembly_id}/{sample_id}.{library_id}.bam",
        reference=METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        tsv=COVERM_BINNING / "{assembly_id}/{sample_id}.{library_id}_genome.tsv",
        euk=COVERM_BINNING / "{assembly_id}/{sample_id}.{library_id}_genome_euk.tsv",
    conda:
        "../envs/binning.yml"
    log:
        COVERM_BINNING / "{assembly_id}/{sample_id}.{library_id}.log",
    params:
        methods=params["binning"]["coverm"]["genome"]["methods"],
        min_covered_fraction=params["binning"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["binning"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} 2> {log}

        coverm genome \
            --separator {params.separator} \
            --bam-files {input.bam} \
            --methods relative_abundance count mean covered_fraction \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.euk} 2>> {log}
        """


rule binning_coverm_aggregate:
    input:
        tsvs=[
            COVERM_BINNING / f"{assembly_id}/{sample_id}.{library_id}_genome.tsv"
            for assembly_id in ["{assemby_id}"]
            for sample_id, library_id in (
                samples[samples.assembly_id == assembly_id][
                    ["sample_id", "library_id"]
                ].values.tolist()
            )
        ],
    output:
        tsv=COVERM_BINNING / "{assembly_id}.coverm_genome.tsv",
    log:
        COVERM_BINNING / "{assembly_id}.coverm_genome.log",
    conda:
        "../envs/binning.yml"
    params:
        input_dir=lambda wildcards: COVERM_BINNING / f"{wildcards.assembly_id}",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule binning_coverm_aggregate_all:
    input:
        tsvs=[
            COVERM_BINNING / f"{assembly_id}.coverm_genome.tsv"
            for assembly_id in samples.assembly_id
        ],
