rule preprocess__bracken__recompute:
    input:
        database=lambda w: features["databases"]["kraken2"][w.kraken2_db],
        report=PRE_KRAKEN2 / "{kraken2_db}" / "{sample_id}.k2report",
    output:
        bracken=touch(
            PRE_BRACKEN / "{kraken2_db}" / "recompute" / "{sample_id}.{level}.bracken"
        ),
        report=touch(
            PRE_BRACKEN / "{kraken2_db}" / "recompute" / "{sample_id}.{level}.report"
        ),
    log:
        PRE_BRACKEN / "{kraken2_db}" / "recompute" / "{sample_id}.{level}.log",
    conda:
        "../../environments/bracken.yml"
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
            -w {output.report} \
            -l {params.level} \
            {params.extra} \
        2> {log} 1>&2
        """


rule preprocess__bracken__combine:
    """Combine all the bracken outputs for a single database"""
    input:
        lambda w: [
            PRE_BRACKEN / w.kraken2_db / "recompute" / f"{sample_id}.{w.level}.bracken"
            for sample_id in SAMPLES
        ],
    output:
        PRE_BRACKEN / "{kraken2_db}" / "combine" / "{level}.tsv",
    log:
        PRE_BRACKEN / "{kraken2_db}" / "combine" / "{level}.log",
    conda:
        "../../environments/bracken.yml"
    shell:
        """
        combine_bracken_outputs.py \
            --files {input} \
            --output {output} \
        2> {log} 1>&2
        """


rule preprocess__bracken__alpha_diversity:
    """Calculate alpha diversity metrics for all samples, metrics and levels at once"""
    input:
        lambda w: [
            PRE_BRACKEN / w.kraken2_db / "recompute" / f"{sample_id}.{level}.bracken"
            for sample_id in SAMPLES
            for level in "DPCOFGS"
        ],
    output:
        PRE_BRACKEN / "{kraken2_db}" / "alpha.tsv",
    log:
        PRE_BRACKEN / "{kraken2_db}" / "alpha.log",
    conda:
        "../../environments/bracken.yml"
    threads: 8
    shell:
        """
        ( parallel \
            --tag \
            --keep-order \
            --jobs {threads} \
            alpha_diversity.py \
                --filename {{1}} \
                --alpha {{2}} \
        ::: {input} \
        ::: Sh BP Si ISi F \
        | sed 's/: /\\t/' \
        ) > {output} \
        2> {log}
        """


rule preprocess__bracken__beta_diversity:
    """Compute beta diversity metrics for all samples and one level"""
    input:
        lambda w: [
            PRE_BRACKEN / w.kraken2_db / "recompute" / f"{sample_id}.{w.level}.bracken"
            for sample_id in SAMPLES
        ],
    output:
        PRE_BRACKEN / "{kraken2_db}" / "beta.{level}.tsv",
    log:
        PRE_BRACKEN / "{kraken2_db}" / "beta.{level}.log",
    conda:
        "../../environments/bracken.yml"
    shell:
        """
        beta_diversity.py \
            --input-files {input} \
            --type bracken \
        > {output} \
        2> {log}
        """


rule preprocess__bracken__all:
    """Get the combined bracken results for all databases"""
    input:
        [
            PRE_BRACKEN / kraken2_db / "recompute" / f"{sample_id}.{level}.bracken"
            for kraken2_db in KRAKEN2_DBS
            for sample_id in SAMPLES
            for level in "DPCOFGS"
        ],
        [
            PRE_BRACKEN / kraken2_db / "combine" / f"{level}.tsv"
            for kraken2_db in KRAKEN2_DBS
            for level in "DPCOFGS"
        ],
        [PRE_BRACKEN / kraken2_db / "alpha.tsv" for kraken2_db in KRAKEN2_DBS],
        [
            PRE_BRACKEN / kraken2_db / f"beta.{level}.tsv"
            for kraken2_db in KRAKEN2_DBS
            for level in "DPCOFGS"
        ],
