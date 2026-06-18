# Input prep: MVP needs one assembly + one read-set per sample ----


rule viruses__mvp__decompress_assembly:
    """Decompress a sample's single-sample assembly: MVP can't parse gzip-ed FASTA"""
    input:
        ASMB_MEGAHIT / "{sample_id}.fa.gz",
    output:
        VIR_MVP_INPUT / "{sample_id}.fa",
    log:
        VIR_MVP_INPUT / "{sample_id}.decompress.log",
    shell:
        """
        gzip --decompress --stdout {input} > {output} 2> {log}
        """


use rule concatenate__gzip_text_files as viruses__mvp__concatenate_forwards with:
    input:
        lambda w: [
            PRE_CLEAN / f"{sample_id}.{library_id}_1.fq.gz"
            for sample_id, library_id in get_libraries_from_sample(w)
        ],
    output:
        VIR_MVP_INPUT / "{sample_id}_R1.fastq.gz",
    log:
        VIR_MVP_INPUT / "{sample_id}_R1.log",


use rule concatenate__gzip_text_files as viruses__mvp__concatenate_reverses with:
    input:
        lambda w: [
            PRE_CLEAN / f"{sample_id}.{library_id}_2.fq.gz"
            for sample_id, library_id in get_libraries_from_sample(w)
        ],
    output:
        VIR_MVP_INPUT / "{sample_id}_R2.fastq.gz",
    log:
        VIR_MVP_INPUT / "{sample_id}_R2.log",


rule viruses__mvp__prepare_metadata:
    """Generate the MVP metadata.tsv: one row per sample, Sample_number = position in SAMPLES"""
    input:
        assemblies=[VIR_MVP_INPUT / f"{sample_id}.fa" for sample_id in SAMPLES],
        forwards=[VIR_MVP_INPUT / f"{sample_id}_R1.fastq.gz" for sample_id in SAMPLES],
        reverses=[VIR_MVP_INPUT / f"{sample_id}_R2.fastq.gz" for sample_id in SAMPLES],
    output:
        VIR_MVP / "metadata.tsv",
    log:
        VIR_MVP / "metadata.log",
    params:
        sample_ids=SAMPLES,
        input_dir=VIR_MVP_INPUT,
    shell:
        """
        ( echo -e "Sample_number\tSample\tAssembly_Path\tRead_Path"
        i=0
        for sample_id in {params.sample_ids} ; do
            i=$((i + 1))
            echo -e "${{i}}\t${{sample_id}}\t{params.input_dir}/${{sample_id}}.fa\t{params.input_dir}/${{sample_id}}_R1.fastq.gz"
        done ) > {output} 2> {log}
        """
