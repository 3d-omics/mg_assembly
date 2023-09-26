include: "functions.smk"


rule bin_metawrap_prepare_one:
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
        "metawrap.yml"
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


rule bin_metawrap_prepare:
    input:
        [
            BIN_METAWRAP / f"{assembly_id}/work_files/{assembly_id}.bam"
            for assembly_id in samples.assembly_id
        ],


rule bin_metawrap_binning_one:
    """Run metawrap over one assembly group
    Note: metawrap works with fastq files, but we can trick it into working by
    creating mock fastq and reference files. h/t: Raphael Eisenhofer
    Note2: metawrap is rotten. It is written in py27 and has a lot unmetable dependencies in conda.
    Using a docker container instead.
    """
    input:
        bam=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}.bam",
        forward_=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}_1.fastq",
        reverse_=BIN_METAWRAP / "{assembly_id}/work_files/{assembly_id}_2.fastq",
        assembly=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        metabat2_bins=directory(BIN_METAWRAP / "{assembly_id}/metabat2_bins"),
        maxbin2_bins=directory(BIN_METAWRAP / "{assembly_id}/maxbin2_bins"),
        concoct_bins=directory(BIN_METAWRAP / "{assembly_id}/concoct_bins"),
    log:
        BIN_METAWRAP / "{assembly_id}.log",
    container:
        "docker://quay.io/biocontainers/metawrap-mg:1.3.0--hdfd78af_1"
    threads: 8
    params:
        min_length=params["bin"]["metawrap"]["binning"]["min_length"],
        extra=params["bin"]["metawrap"]["binning"]["extra"],
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


rule bin_metawrap:
    input:
        [
            BIN_METAWRAP / f"{assembly_id}/{binner}"
            for assembly_id in samples.assembly_id
            for binner in ["concoct_bins", "maxbin2_bins", "metabat2_bins"]
        ],
