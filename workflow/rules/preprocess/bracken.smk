include: "bracken_functions.smk"


rule preprocess__bracken__recompute:
    """Recompute kraken2 reports and counts and for all samples and levels"""
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
        ENVS / "bracken.yml"
    params:
        read_length=params["preprocess"]["bracken"]["read_length"],
        threshold=params["preprocess"]["bracken"]["threshold"],
        level=lambda w: w.level,
    shell:
        """
        if [ ! -s {input.report} ] ; then
            echo "Empty kraken2 report. Skipping" 2> {log} 1>&2
            echo -e "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads" > {output.bracken}
            exit 0
        fi

        {{
            bracken \
                -d {input.database} \
                -i {input.report} \
                -o {output.bracken} \
                -w {output.report} \
                -l {params.level} \
                -r {params.read_length} \
                -t {params.threshold} \
            2> {log} 1>&2
        }} || {{
            echo -e "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads" > {output.bracken}
        }}
        """


rule preprocess__bracken__recompute__all:
    """Run preprocess__bracken__recompute for all databases, samples and levels"""
    input:
        [
            PRE_BRACKEN / kraken2_db / "recompute" / f"{sample_id}.{level}.bracken"
            for kraken2_db in KRAKEN2_DBS
            for sample_id in SAMPLES
            for level in ALL_TAXONOMY_LEVELS
        ],


rule preprocess__bracken__report:
    """Move a species report to a different folder for MultiQC"""
    input:
        PRE_BRACKEN / "{kraken2_db}" / "recompute" / "{sample_id}.S.report",
    output:
        PRE_BRACKEN / "{kraken2_db}" / "report" / "{sample_id}.report",
    log:
        PRE_BRACKEN / "{kraken2_db}" / "report" / "{sample_id}.log",
    conda:
        ENVS / "bracken.yml"
    shell:
        """
        cp --verbose {input} {output} 2> {log} 1>&2
        """


rule preprocess__bracken__report__all:
    """Run preprocess__bracken__report for all databases and samples"""
    input:
        [
            PRE_BRACKEN / kraken2_db / "report" / f"{sample_id}.report"
            for kraken2_db in KRAKEN2_DBS
            for sample_id in SAMPLES
        ],


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
        ENVS / "bracken.yml"
    shell:
        """
        combine_bracken_outputs.py \
            --files {input} \
            --output {output} \
        2> {log} 1>&2
        """


rule preprocess__bracken__combine__all:
    """Run preprocess__bracken__combine for all databases and levels"""
    input:
        [
            PRE_BRACKEN / kraken2_db / "combine" / f"{level}.tsv"
            for kraken2_db in KRAKEN2_DBS
            for level in ALL_TAXONOMY_LEVELS
        ],


rule preprocess__bracken__alpha_diversity:
    """Calculate alpha diversity metrics for all samples, metrics and levels at once"""
    input:
        lambda w: [
            PRE_BRACKEN / w.kraken2_db / "recompute" / f"{sample_id}.{w.level}.bracken"
            for sample_id in SAMPLES
        ],
    output:
        PRE_BRACKEN / "{kraken2_db}" / "alpha.{level}.tsv",
    log:
        PRE_BRACKEN / "{kraken2_db}" / "alpha.{level}.log",
    conda:
        ENVS / "bracken.yml"
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
        | grep -v loading \
        ) > {output} \
        2> {log}
        """


rule preprocess__bracken__alpha_diversity__all:
    """Run preprocess__bracken__alpha_diversity for all databases and levels"""
    input:
        [
            PRE_BRACKEN / kraken2_db / f"alpha.{level}.tsv"
            for kraken2_db in KRAKEN2_DBS
            for level in ALL_TAXONOMY_LEVELS
        ],


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
        ENVS / "bracken.yml"
    shell:
        """
        beta_diversity.py \
            --input-files {input} \
            --type bracken \
        > {output} \
        2> {log}
        """


rule preprocess__bracken__beta_diversity__all:
    """Run preprocess__bracken__beta_diversity for all databases and some levels

    Note: Beta diversity is only computed for Order, Family, Genus and Species
    """
    input:
        [
            PRE_BRACKEN / kraken2_db / f"beta.{level}.tsv"
            for kraken2_db in KRAKEN2_DBS
            for level in LOW_TAXONOMY_LEVELS
        ],


rule preprocess__bracken__all:
    """Get the combined bracken results for all databases"""
    input:
        rules.preprocess__bracken__recompute__all.input,
        rules.preprocess__bracken__report__all.input,
        rules.preprocess__bracken__combine__all.input,
        rules.preprocess__bracken__alpha_diversity__all.input,
        rules.preprocess__bracken__beta_diversity__all.input,
