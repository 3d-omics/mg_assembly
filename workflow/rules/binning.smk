rule binning_prepare_one:
    input:
        bams=get_bams_for_metawrap_binning_prepare,
    output:
        forward_=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}_1.fastq",
        reverse_=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}_2.fastq",
        bwt_index=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}.fa.bwt",
        bam=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}.bam",
    log:
        METAWRAP_BINNING / "{assembly_id}.prepare.log",
    conda:
        "../envs/binning.yml"
    params:
        n=get_number_of_libraries_in_binning,
    shell:
        """
        echo "@" > {output.forward_} 2> {log}
        echo "@" > {output.reverse_} 2>> {log}
        touch {output.bwt_index} 2>> {log} 1>&2
        if [[ {params.n} -eq 1 ]] ; then
            samtools view -F 4 -1 {input.bams} > {output.bam} 2>> {log}
        else
            (samtools merge \
                -u \
                -o /dev/stdout \
                {input.bams} \
            | samtools view \
                -o {output.bam} \
                -F 4 \
                -1 \
            ) 2>> {log}
        fi
        """


rule binning_prepare:
    input:
        [
            METAWRAP_BINNING / f"{assembly_id}/work_files/{assembly_id}.bam"
            for assembly_id in samples.assembly_id
        ],


rule binning_binning_one:
    """Run metawrap over one assembly group
    Note: metawrap works with fastq files, but we can trick it into working by
    creating mock fastq and reference files. h/t: Raphael Eisenhofer
    Note2: metawrap is rotten. It is written in py27 and has a lot unmetabel dependencies in conda.
    Using a singularity container instead.
    """
    input:
        bam=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}.bam",
        forward_=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}_1.fastq",
        reverse_=METAWRAP_BINNING / "{assembly_id}/work_files/{assembly_id}_2.fastq",
        assembly=MEGAHIT_RENAMING / "{assembly_id}.fa",
    output:
        metabat2_bins=directory(METAWRAP_BINNING / "{assembly_id}/metabat2_bins"),
        maxbin2_bins=directory(METAWRAP_BINNING / "{assembly_id}/maxbin2_bins"),
        concoct_bins=directory(METAWRAP_BINNING / "{assembly_id}/concoct_bins"),
    log:
        METAWRAP_BINNING / "{assembly_id}.log",
    singularity:
        "https://depot.galaxyproject.org/singularity/metawrap-mg:1.3.0--hdfd78af_1"
    threads: 8
    params:
        min_length=params["binning"]["metawrap_binning"]["min_length"],
        extra=params["binning"]["metawrap_binning"]["extra"],
        out_folder=lambda wildcards: METAWRAP_BINNING / f"{wildcards.assembly_id}",
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        metawrap binning \
            -o {params.out_folder} \
            -t {threads} \
            -m $(({resources.mem_mb} / 1024)) \
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


rule binning_binning:
    input:
        [
            METAWRAP_BINNING / f"{assembly_id}/{binner}"
            for assembly_id in samples.assembly_id
            for binner in ["concoct_bins", "maxbin2_bins", "metabat2_bins"]
        ],


rule binning_refinement_one:
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
    threads: 16
    params:
        completeness=params["binning"]["metawrap_bin_refinement"]["completeness"],
        contamination=params["binning"]["metawrap_bin_refinement"]["contamination"],
        extra=params["binning"]["metawrap_bin_refinement"]["extra"],
        output_prefix=compose_metawrap_working_folder,
    resources:
        mem_mb=8 * 1024,
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


rule binning_refinement_all:
    input:
        [
            METAWRAP_REFINEMENT / f"{assembly_id}.contigs"
            for assembly_id in samples.assembly_id
        ],


rule binning_renaming_one:
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
            bin_name=$(basename $bin .fa); \
            seqtk rename $bin {params.assembly_id}.${{bin_name}}. ; \
        done > {output.fa}) 2> {log}
        """


rule binning_renaming_all:
    input:
        [METAWRAP_RENAMING / f"{assembly_id}.fa" for assembly_id in samples.assembly_id],


rule binning_index_one:
    input:
        bins=METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        mock=touch(BOWTIE2_INDEXES_BINNING / "{assembly_id}"),
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
            {input.bins} \
            {output.mock} \
        2> {log} 1>&2
        """


rule binning_bowtie2_one:
    input:
        mock=BOWTIE2_INDEXES_BINNING / "{assembly_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        reference=METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        cram=BOWTIE2_BINNING / "{assembly_id}.{sample_id}.{library_id}.cram",
    log:
        BOWTIE2_BINNING / "{assembly_id}.{sample_id}.{library_id}.log",
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


rule binning_bowtie2:
    input:
        [
            BOWTIE2_BINNING / f"{assembly_id}.{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule binning_cram_to_bam_one:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=BOWTIE2_BINNING / "{assembly_id}.{sample_id}.{library_id}.cram",
        crai=BOWTIE2_BINNING / "{assembly_id}.{sample_id}.{library_id}.cram.crai",
        reference=METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        bam=temp(BOWTIE2_BINNING / "{assembly_id}.{sample_id}.{library_id}.bam"),
    log:
        BOWTIE2_BINNING / "{assembly_id},{sample_id}.{library_id}.bam.log",
    conda:
        "../envs/binning.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=4 * 1024,
    shell:
        """
        samtools view \
            -F 4 \
            --threads {threads} \
            --reference {input.reference} \
            --output {output.bam} \
            --fast \
            {input.cram} \
        2> {log}
        """


rule binning_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=BOWTIE2_BINNING / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        tsv=COVERM_BINNING / "genome/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/binning.yml"
    log:
        COVERM_BINNING / "genome/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["assembly"]["coverm"]["genome"]["methods"],
        min_covered_fraction=params["assembly"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["assembly"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} 2> {log}
        """


rule binning_coverm_genome:
    input:
        tsvs=[
            COVERM_BINNING / f"genome/{assembly_id}.{sample_id}.{library_id}.tsv"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
    output:
        tsv=COVERM_BINNING / "genome.tsv",
    log:
        COVERM_BINNING / "genome.log",
    conda:
        "../envs/assembly.yml"
    params:
        input_dir=COVERM_BINNING / "genome/",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule binning_coverm_contig_one:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=BOWTIE2_BINNING / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        tsv=COVERM_BINNING / "contig/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/binning.yml"
    log:
        COVERM_BINNING / "contig/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["binning"]["coverm"]["contig"]["methods"],
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output.tsv} 2> {log}
        """


rule binning_coverm_contig:
    input:
        tsvs=[
            COVERM_BINNING / f"contig/{assembly_id}.{sample_id}.{library_id}.tsv"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
    output:
        tsv=COVERM_BINNING / "contig.tsv",
    log:
        COVERM_BINNING / "contig.log",
    conda:
        "../envs/binning.yml"
    params:
        input_dir=COVERM_BINNING / f"contig",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule binning_quast_one:
    """Run quast over one assembly group"""
    input:
        METAWRAP_RENAMING / "{assembly_id}.fa",
    output:
        directory(BINNING_QUAST / "{assembly_id}"),
    log:
        BINNING_QUAST / "{assembly_id}.log",
    conda:
        "../envs/binning.yml"
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


rule binning_quast_all:
    """Run quast over all assembly groups"""
    input:
        [BINNING_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule binning:
    input:
        rules.binning_coverm_contig.output,
        rules.binning_coverm_genome.output,
        rules.binning_quast_all.input,
