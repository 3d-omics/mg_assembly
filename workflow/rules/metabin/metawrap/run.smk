rule metawrap_refinement_one:
    input:
        metabat2_bins=BIN_METAWRAP / "{assembly_id}/metabat2_bins",
        maxbin2_bins=BIN_METAWRAP / "{assembly_id}/maxbin2_bins",
        concoct_bins=BIN_METAWRAP / "{assembly_id}/concoct_bins",
    output:
        stats=METABIN_METAWRAP / "{assembly_id}_bins.stats",
        contigs=METABIN_METAWRAP / "{assembly_id}_bins.contigs",
        working_folder=directory(METABIN_METAWRAP / "{assembly_id}"),
        bins_folder=directory(METABIN_METAWRAP / "{assembly_id}_bins"),
    log:
        METABIN_METAWRAP / "{assembly_id}.log",
    container:
        "docker://quay.io/biocontainers/metawrap-mg:1.3.0--hdfd78af_1"
    threads: 16
    params:
        completeness=params["metabin"]["metawrap"]["refinement"]["completeness"],
        contamination=params["metabin"]["metawrap"]["refinement"]["contamination"],
        extra=params["metabin"]["metawrap"]["refinement"]["extra"],
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


rule metawrap_renaming_one:
    """
    Note: doing this separatedly from the binning step because we need seqtk and it is outside the metawrap singularity container
    """
    input:
        bin_folder=METABIN_METAWRAP / "{assembly_id}_bins",
    output:
        fa=METAWRAP_RENAME / "{assembly_id}.fa",
    log:
        ASSEMBLY_RENAME / "{assembly_id}.log",
    conda:
        "../envs/metabin.yml"
    params:
        assembly_id=lambda wildcards: f"{wildcards.assembly_id}",
    shell:
        """
        (for bin in {input.bin_folder}/*.fa ; do
            bin_name=$(basename $bin .fa); \
            seqtk rename $bin {params.assembly_id}.${{bin_name}}. ; \
        done > {output.fa}) 2> {log}
        """


rule metawrap_renaming_all:
    input:
        [ASSEMBLY_RENAME / f"{assembly_id}.fa" for assembly_id in samples.assembly_id],


rule metawrap:
    input:
        [
            METABIN_METAWRAP / f"{assembly_id}_bins.contigs"
            for assembly_id in samples.assembly_id
        ],
