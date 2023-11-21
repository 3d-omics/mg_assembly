include: "_functions.smk"
include: "fastp.smk"
include: "bowtie2.smk"
include: "fastqc.smk"
include: "nonpareil.smk"
include: "singlem.smk"
include: "coverm.smk"
include: "kraken2.smk"
include: "samtools.smk"


rule preprocess__eval:
    """Run the evaluation of the preprocessing steps."""
    input:
        rules.preprocess__fastqc.input,
        rules.preprocess__coverm.output,
        rules.preprocess__samtools.input,
        rules.preprocess__kraken2.output,


rule prerocess__eval_with_singlem:
    """Run the evaluation of the preprocessing steps + SingleM."""
    input:
        rules.preprocess__eval.input,
        rules.preprocess__singlem.output,


rule preprocess__eval_with_nonpareil:
    """Run the evaluation of the preprocessing steps + Nonpareil."""
    input:
        rules.preprocess__eval.input,
        rules.preprocess__nonpareil.output,


rule preprocess__run:
    """Run the preprocessing steps, without the evaluation ones"""
    input:
        rules.preprocess__bowtie2__extract_nonhost.input,


rule preprocess:
    """Run the preprocessing steps, included the evaluation ones"""
    input:
        rules.preprocess__run.input,
        rules.preprocess__eval.input,
