rule reads_eval_fastqc:
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.zip"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule reads_eval:
    input:
        rules.reads_eval_fastqc.input,
