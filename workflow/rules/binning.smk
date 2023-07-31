include: "binning/vamb.smk"
# include: "binning/concoct.smk"
include: "binning/metabat2.smk"


# rule binning_metabat2_one:
#     input:
#         assembly=MEGAHIT_RENAMING / "{assembly_id}.fa",
#     output:

#     shell:


# rule binning_maxbin2_one:
#     input:

#     output:

#     shell:


# rule binning_magscot_prodigal_one:
#     input:
#         assembly=MEGAHIT_RENAMING / "{assembly_id}.fa",
#     output:
#         proteins = MAGSCOT_PRODIGAL / "{assembly_id}.faa",
#         genes = MAGSCOT_PRODIGAL / "{assembly_id}.ffn",
#         txt = MAGSCOT_PRODIGAL / "{assembly_id}.txt",
#     log:
#         MAGSCOT_PRODIGAL / "{assembly_id}.log",
#     shell:
#         """
#         prodigal \
#             -i {input.assembly} \
#             -o {output.txt} \
#             -d {output.genes}
#             -a {output.proteins} \
#             -p meta \
#         2> {log} 1>&2
#         """


# rule binning_magscot_hmmsearch_pfam_one:
#     input:
#         proteins = MAGSCOT_PRODIGAL / "{assembly_id}.faa",
#         hmm = features["binning"]["magscot"]["pfam_hmm"],
#     output:
#         out = MAGSCOT_HMM / "{assembly_id}.pfam.out",
#         tblout = MAGSCOT_HMM / "{assembly_id}.pfam.tblout",
#     log:
#         MAGSCOT_HMM / "{assembly_id}.pfam.log",
#     conda:
#         "../envs/binning.yml"
#     threads:
#         16
#     shell:
#         """
#         hmmsearch \
#             -o example.hmm.tigr.out \
#             --tblout example.hmm.tigr.hit.out \
#             --noali \
#             --notextw \
#             --cut_nc \
#             --cpu {threads} \
#             $MAGScoT_folder/hmm/gtdbtk_rel207_tigrfam.hmm \
#             {input.proteins} \
#         2> {log} 1>&2
#         """


# rule binning_magscot_hmmsearch_tigr_one:
#     input:
#         proteins = MAGSCOT_PRODIGAL / "{assembly_id}.faa",
#         hmm = features["binning"]["magscot"]["tigr"],
#     output:
#         out = MAGSCOT_HMM / "{assembly_id}.tigr.out",
#         tblout = MAGSCOT_HMM / "{assembly_id}.tigr.tblout",
#     log:
#         MAGSCOT_HMM / "{assembly_id}.tigr.log",
#     conda:
#         "../envs/binning.yml"
#     threads:
#         16
#     shell:
#         """
#         hmmsearch \
#             -o example.hmm.tigr.out \
#             --tblout example.hmm.tigr.hit.out \
#             --noali \
#             --notextw \
#             --cut_nc \
#             --cpu {threads} \
#             $MAGScoT_folder/hmm/gtdbtk_rel207_tigrfam.hmm \
#             {input.proteins} \
#         2> {log} 1>&2
#         """


# rule binning_magscot_join_hmm_one:
#     input:
#         tigr_out = MAGSCOT_HMM / "{assembly_id}.tigr.out",
#         pfam_out = MAGSCOT_HMM / "{assembly_id}.pfam.out",
#     output:
#         merged = MAGSCOT_HMM / "{assembly_id}.out",
#     log:
#         MAGSCOT_HMM / "{assembly_id}.log",
#     conda:
#         "../envs/binning.yml"
#     shell:
#         """
#         grep -v "^#" {input.tigr_out} | awk '{print $1"\t"$3"\t"$5}' >  {output.merged} 2>  {log}
#         grep -v "^#" {input.pfam_out} | awk '{print $1"\t"$4"\t"$5}' >> {output.merged} 2>> {log}
#         """


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
