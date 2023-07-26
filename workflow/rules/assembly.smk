rule assembly_megahit_one:
    """Run megahit over one sample, merging all libraries in the process

    Note: the initial rm -rf is to delete the folder that snakemake creates.
    megahit refuses to overwrite an existing folder
    """
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        assembly_folder=directory(MEGAHIT / "{assembly_id}"),
        assembly_fasta=MEGAHIT / "{assembly_id}/final.contigs.fa",
    log:
        log=MEGAHIT / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 16
    resources:
        mem_mb=64 * 1024,
    params:
        min_contig_len=params["assembly"]["megahit"]["min_contig_len"],
        extra=params["assembly"]["megahit"]["extra"],
        forwards=aggregate_forwards_for_megahit,
        reverses=aggregate_reverses_for_megahit,
    shell:
        """
        megahit \
            --num-cpu-threads {threads} \
            --min-contig-len {params.min_contig_len} \
            --verbose \
            --force \
            -1 {params.forwards} \
            -2 {params.reverses} \
            {params.extra} \
            --out-dir {output.assembly_folder} \
        2> {log} 1>&2
        """


rule assembly_megahit_all:
    """Run megahit over all groups"""
    input:
        [MEGAHIT / f"{assembly_id}.contigs.fa" for assembly_id in samples.assembly_id],


rule assembly_quast_one:
    """Run quast over one assembly group"""
    input:
        MEGAHIT / "{assembly_id}/final.contigs.fa",
    output:
        directory(QUAST / "{assembly_id}"),
    log:
        QUAST / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 4
    params:
        extra=params["assembly"]["quast"]["extra"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """


rule assembly_quast_all:
    """Run quast over all assembly groups"""
    input:
        [QUAST / f"{assembly_id}" for assembly_id in samples.assembly_id],


rule assembly:
    """Run all the assemblies"""
    input:
        rules.assembly_quast_all.input,


rule assembly_bowtie2_build_one:
    """
    Index megahit assembly
    """
    input:
        contigs=MEGAHIT / "{assembly_id}" / "final.contigs.fa",
    output:
        mock=touch(BOWTIE2_INDEXES_ASSEMBLY / "{assembly_id}"),
    log:
        BOWTIE2_INDEXES_ASSEMBLY / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
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


rule assembly_bowtie2_build_all:
    input:
        [
            BOWTIE2_INDEXES_ASSEMBLY / f"{assembly_id}"
            for assembly_id in samples.assembly_id
        ],


rule assembly_bowtie2_one:
    input:
        mock=BOWTIE2_INDEXES_ASSEMBLY / "{assembly_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        reference=MEGAHIT / "{assembly_id}" / "final.contigs.fa",
    output:
        cram=BOWTIE2_ASSEMBLY / "{assembly_id}/{sample_id}.{library_id}.cram",
    log:
        log=BOWTIE2_ASSEMBLY / "{assembly_id}/{sample_id}.{library_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 24
    params:
        extra=params["assembly"]["bowtie2"]["extra"],
        samtools_mem=params["assembly"]["samtools"]["mem"],
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


rule assembly_bowtie2_all:
    input:
        [
            BOWTIE2_ASSEMBLY / f"{assembly_id}/{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in samples[
                ["assembly_id", "sample_id", "library_id"]
            ].values.tolist()
        ],


rule assembly_merge_bams_one:
    input:
        crams=get_crams_to_merge,
        reference=MEGAHIT / "{assembly_id}" / "final.contigs.fa",
    output:
        bam=BOWTIE2_ASSEMBLY / "{assembly_id}.bam",
    log:
        log=BOWTIE2_ASSEMBLY / "{assembly_id}.log",
    conda:
        "../envs/assembly.yml"
    threads: 24
    params:
        n=get_number_of_libraries_in_assembly,
    shell:
        """
        if [ {params.n} -eq 1 ] ; then
            samtools view \
                --bam \
                --reference {input.reference} \
                {input.crams} \
            > {output.bam} \
            2> {log}
        else
            samtools merge \
                -@ {threads} \
                -l 1 \
                -o {output.bam} \
                {input.crams} \
            2> {log} 1>&2
        fi
        """


rule assembly_merge_bams_all:
    input:
        [BOWTIE2_ASSEMBLY / f"{assembly_id}.bam" for assembly_id in samples.assembly_id],


rule assembly_metawrap_binning_prepare_one:
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
        "../envs/assembly.yml"
    shell:
        """
        echo "@" > {output.forward_}
        echo "@" > {output.reverse_}
        touch {output.bwt_index}
        ln {input.bam} {output.bam}
        """


rule assembly_metawrap_binning_prepare_all:
    input:
        [METAWRAP_BINNING / f"{assembly_id}" for assembly_id in samples.assembly_id],


rule assembly_metawrap_binning_one:
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
        min_length=params["assembly"]["metawrap_binning"]["min_length"],
        extra=params["assembly"]["metawrap_binning"]["extra"],
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


rule assembly_metawrap_binning_all:
    input:
        [
            METAWRAP_BINNING / f"{assembly_id}/{binner}"
            for assembly_id in samples.assembly_id
            for binner in ["concoct_bins", "maxbin2_bins", "metabat2_bins"]
        ],


rule assembly_metawrap_bin_refinement_one:
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
        completeness=params["assembly"]["metawrap_bin_refinement"]["completeness"],
        contamination=params["assembly"]["metawrap_bin_refinement"]["contamination"],
        extra=params["assembly"]["metawrap_bin_refinement"]["extra"],
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
            `# {METAWRAP_REFINEMENT}/{wildcards.assembly_id}/metawrap_{params.completeness}_{params.contamination}_bins.stats` \
            {output.stats} \
        2>> {log} 1>&2

        cp \
            {params.output_prefix}.contigs \
            `# {METAWRAP_REFINEMENT}/{wildcards.assembly_id}/metawrap_{params.completeness}_{params.contamination}_bins.contigs` \
            {output.contigs} \
        2>> {log} 1>&2

        cp --recursive \
            {params.output_prefix} \
            {output.bins_folder} \
        2>> {log} 1>&2
        """


rule assembly_metawrap_bin_refinement_all:
    input:
        [
            METAWRAP_REFINEMENT / f"{assembly_id}.contigs"
            for assembly_id in samples.assembly_id
        ],


rule assembly_metawrap_renaming_one:
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
        "../envs/assembly.yml"
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


rule assembly_metawrap_renaming_all:
    input:
        [METAWRAP_RENAMING / f"{assembly_id}.fa" for assembly_id in samples.assembly_id],
