include: "metawrap_functions.smk"


rule metawrap_prepare_one:
    input:
        bams=get_bams_for_metawrap_metawrap_prepare,
    output:
        forward_=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}_1.fastq",
        reverse_=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}_2.fastq",
        bwt_index=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}.fa.bwt",
        bam=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}.bam",
    log:
        BIN_METAWRAP / "{assembly_id}.prepare.log",
    conda:
        "../../envs/binning/metawrap.yml"
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


rule metawrap_prepare:
    input:
        [
            BIN_METAWRAP / f"{assembly_id}/work_files/{assembly_id}.bam"
            for assembly_id in samples.assembly_id
        ],


rule metawrap_binning_one:
    """Run metawrap over one assembly group
    Note: metawrap works with fastq files, but we can trick it into working by
    creating mock fastq and reference files. h/t: Raphael Eisenhofer
    Note2: metawrap is rotten. It is written in py27 and has a lot unmetabel dependencies in conda.
    Using a singularity container instead.
    """
    input:
        bam=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}.bam",
        forward_=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}_1.fastq",
        reverse_=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}_2.fastq",
        assembly=ASSEMBLY_RENAME / "{assembly_id}.fa",
    output:
        metabat2_bins=directory(BIN_METAWRAP / "{assembly_id}/metabat2_bins"),
        maxbin2_bins=directory(BIN_METAWRAP / "{assembly_id}/maxbin2_bins"),
        concoct_bins=directory(BIN_METAWRAP / "{assembly_id}/concoct_bins"),
    log:
        BIN_METAWRAP / "{assembly_id}.log",
    singularity:
        "https://depot.galaxyproject.org/singularity/metawrap-mg:1.3.0--hdfd78af_1"
    threads: 8
    params:
        min_length=params["binning"]["metawrap"]["binning"]["min_length"],
        extra=params["binning"]["metawrap"]["binning"]["extra"],
        out_folder=lambda wildcards: BIN_METAWRAP / f"{wildcards.assembly_id}",
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


rule metawrap_binning:
    input:
        [
            METABIN_METAWRAP / f"{assembly_id}/{binner}"
            for assembly_id in samples.assembly_id
            for binner in ["concoct_bins", "maxbin2_bins", "metabat2_bins"]
        ],


# rule metawrap_refinement_one:
#     input:
#         metabat2_bins=BIN_METAWRAP / "{assembly_id}/metabat2_bins",
#         maxbin2_bins=BIN_METAWRAP / "{assembly_id}/maxbin2_bins",
#         concoct_bins=BIN_METAWRAP / "{assembly_id}/concoct_bins",
#     output:
#         stats=METABIN_METAWRAP / "{assembly_id}_bins.stats",
#         contigs=METABIN_METAWRAP / "{assembly_id}_bins.contigs",
#         working_folder=directory(METABIN_METAWRAP / "{assembly_id}"),
#         bins_folder=directory(METABIN_METAWRAP / "{assembly_id}_bins"),
#     log:
#         METABIN_METAWRAP / "{assembly_id}.log",
#     singularity:
#         "https://depot.galaxyproject.org/singularity/metawrap-mg:1.3.0--hdfd78af_1"
#     threads: 16
#     params:
#         completeness=params["binning"]["metawrap"]["refinement"]["completeness"],
#         contamination=params["binning"]["metawrap"]["refinement"]["contamination"],
#         extra=params["binning"]["metawrap"]["refinement"]["extra"],
#         output_prefix=compose_metawrap_working_folder,
#     resources:
#         mem_mb=8 * 1024,
#     shell:
#         """
#         metawrap bin_refinement \
#             -m $(({resources.mem_mb} / 1024)) \
#             -o {output.working_folder} \
#             -t {threads} \
#             -A {input.metabat2_bins} \
#             -B {input.maxbin2_bins} \
#             -C {input.concoct_bins} \
#             -c {params.completeness} \
#             -x {params.contamination} \
#             {params.extra} \
#         2> {log} 1>&2
#         cp \
#             {params.output_prefix}.stats \
#             {output.stats} \
#         2>> {log} 1>&2
#         cp \
#             {params.output_prefix}.contigs \
#             {output.contigs} \
#         2>> {log} 1>&2
#         cp --recursive \
#             {params.output_prefix} \
#             {output.bins_folder} \
#         2>> {log} 1>&2
#         """
# rule metawrap_renaming_one:
#     """
#     Note: doing this separatedly from the binning step because we need seqtk and it is outside the metawrap singularity container
#     """
#     input:
#         bin_folder=METABIN_METAWRAP / "{assembly_id}_bins",
#     output:
#         fa= METAWRAP_RENAME / "{assembly_id}.fa",
#     log:
#         ASSEMBLY_RENAME / "{assembly_id}.log",
#     conda:
#         "../envs/binning.yml"
#     params:
#         assembly_id=lambda wildcards: f"{wildcards.assembly_id}",
#     shell:
#         """
#         (for bin in {input.bin_folder}/*.fa ; do
#             bin_name=$(basename $bin .fa); \
#             seqtk rename $bin {params.assembly_id}.${{bin_name}}. ; \
#         done > {output.fa}) 2> {log}
#         """
# rule metawrap_renaming_all:
#     input:
#         [ASSEMBLY_RENAME / f"{assembly_id}.fa" for assembly_id in samples.assembly_id],
# rule metawrap:
#     input:
#         [
#             METABIN_METAWRAP / f"{assembly_id}_bins.contigs"
#             for assembly_id in samples.assembly_id
#         ],
