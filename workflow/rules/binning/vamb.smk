include: "vamb_functions.smk"


rule binning_vamb_concatenate_one:
    input:
        assembly=MEGAHIT_RENAMING / "{assembly_id}.fa",
    output:
        concatenated=VAMB / "concatenated" / "{assembly_id}.fa.gz",
    log:
        VAMB / "concatenated/{assembly_id}.log",
    conda:
        "../envs/binning.yml"
    shell:
        """
        concatenate.py \
            {output.concatenated} \
            {input.assembly} \
        2> {log} 1>&2
        """


rule binning_vamb_index_one:
    input:
        concatenated=VAMB / "concatenated" / "{assembly_id}.fa.gz",
    output:
        index=touch(VAMB / "indexes" / "{assembly_id}"),
    log:
        VAMB / "indexes/{assembly_id}.log",
    conda:
        "../envs/binning.yml"
    threads: 24
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.concatenated} \
            {output.index} \
        2> {log} 1>&2
        """


rule binning_vamb_map_one:
    input:
        index=VAMB / "indexes" / "{assembly_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
    output:
        bam=VAMB / "bams" / "{assembly_id}.{sample_id}.{library_id}.bam",
    log:
        VAMB / "bams/{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        "../envs/binning.yml"
    threads: 24
    shell:
        """
        (bowtie2 \
            --threads {threads} \
            --no-unal \
            -x {input.index} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
        | samtools view \
            -F 3584 \
            -b \
            --threads {threads} \
        > {output.bam} \
        ) 2> {log} 1>&2
        """


rule binning_vamb_one:
    input:
        concatenated=VAMB / "concatenated" / "{assembly_id}.fa.gz",
        bams=get_vamb_bams_from_assembly_id,
    output:
        folder=directory(VAMB / "bins" / "{assembly_id}"),
    log:
        VAMB / "bins/{assembly_id}.log",
    conda:
        "../envs/binning.yml"
    params:
        extra="",
    threads: 1
    shell:
        """
        vamb \
            --outdir {output.folder} \
            --fasta {input.concatenated} \
            --bamfiles {input.bams} \
            {params.extra} \
            -p {threads} \
        2> {log} 1>&2
        """


rule binning_vamb_all:
    input:
        [VAMB / "bins" / f"{assembly_id}" for assembly_id in ASSEMBLIES],
