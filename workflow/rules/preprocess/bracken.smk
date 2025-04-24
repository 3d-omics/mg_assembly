rule preprocess__bracken__recompute:
    input:
        database=lambda w: features["databases"]["kraken2"][w.kraken2_db],
        report=PRE_KRAKEN2 / "{kraken2_db}" / "{sample_id}.k2report",
    output:
        bracken=touch(
            PRE_BRACKEN / "{kraken2_db}" / "{sample_id}.{level}.bracken"
        ),
    log:
        PRE_BRACKEN / "{kraken2_db}" / "{sample_id}.{level}.log",
    conda:
        "../../environments/kraken2.yml"
    params:
        extra=params["preprocess"]["kraken2"]["bracken"]["extra"],
        level=lambda w: w.level,
    shell:
        """
        if [ ! -s {input.report} ] ; then
            echo "Empty report. Skipping" 2> {log} 1>&2
            exit 0
        fi

        bracken \
            -d {input.database} \
            -i {input.report} \
            -o {output.bracken} \
            -l {params.level} \
            {params.extra} \
        2> {log} 1>&2
        """


rule preprocess__bracken__combine:
    """Combine all the bracken outputs for a single database"""
    input:
        lambda w: [
            PRE_KRAKEN2 / w.kraken2_db / f"{sample_id}.{w.level}.bracken"
            for sample_id in SAMPLES
        ],
    output:
        PRE_KRAKEN2 / "{kraken2_db}.{level}.tsv",
    log:
        PRE_KRAKEN2 / "{kraken2_db}.{level}.log",
    conda:
        "../../environments/kraken2.yml"
    shell:
        """
        combine_bracken_outputs.py \
            --files {input} \
            --output {output} \
        2> {log} 1>&2
        """


rule preprocess__bracken__all:
    """Get the combined bracken results for all databases"""
    input:
        [
            PRE_KRAKEN2 / f"{kraken2_db}.{level}.tsv"
            for kraken2_db in KRAKEN2_DBS
            for level in "DPCOFGS"
        ],
