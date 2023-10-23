include: "functions.smk"
include: "fastp.smk"
include: "bowtie2.smk"
include: "fastqc.smk"
include: "nonpareil.smk"
include: "singlem.smk"
include: "coverm.smk"
include: "kraken2.smk"
include: "samtools.smk"


rule pre_eval:
    """Run the evaluation of the preprocessing steps."""
    input:
        rules.pre_fastp_fastqc.input,
        rules.pre_coverm.output,
        rules.pre_samtools.input,
        rules.pre_nonhost_fastqc.input,
        rules.pre_kraken2.output,


rule pre_eval_with_singlem:
    """Run the evaluation of the preprocessing steps + SingleM."""
    input:
        rules.pre_eval.input,
        rules.pre_singlem.output,


rule pre_eval_with_nonpareil:
    """Run the evaluation of the preprocessing steps + Nonpareil."""
    input:
        rules.pre_eval.input,
        rules.pre_nonpareil.output,


rule pre_run:
    """Run the preprocessing steps, without the evaluation ones"""
    input:
        rules.pre_bowtie2_extract_nonhost_all.input,


rule pre:
    """Run the preprocessing steps, included the evaluation ones"""
    input:
        rules.pre_run.input,
        rules.pre_eval.input,
