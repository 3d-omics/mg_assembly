include: "bakta_functions.smk"


rule prokaryotes__annotate__bakta:
    """Annotate a single MAG with Bakta"""
    input:
        fasta=PROK_MAGS / "{mag_id}.fa",
        db=features["databases"]["bakta"],
    output:
        embl=PROK_ANN / "bakta" / "{mag_id}.embl.gz",
        faa=PROK_ANN / "bakta" / "{mag_id}.faa.gz",
        ffn=PROK_ANN / "bakta" / "{mag_id}.ffn.gz",
        fna=PROK_ANN / "bakta" / "{mag_id}.fna.gz",
        gbff=PROK_ANN / "bakta" / "{mag_id}.gbff.gz",
        gff3=PROK_ANN / "bakta" / "{mag_id}.gff3.gz",
        hypotheticals_faa=PROK_ANN / "bakta" / "{mag_id}.hypotheticals.faa.gz",
        hypotheticals_tsv=PROK_ANN / "bakta" / "{mag_id}.hypotheticals.tsv.gz",
        inference=PROK_ANN / "bakta" / "{mag_id}.inference.tsv.gz",
        json=PROK_ANN / "bakta" / "{mag_id}.json.gz",
        png=PROK_ANN / "bakta" / "{mag_id}.png.gz",
        svg=PROK_ANN / "bakta" / "{mag_id}.svg",
        tsv=PROK_ANN / "bakta" / "{mag_id}.tsv.gz",
    log:
        PROK_ANN / "bakta" / "{mag_id}.log",
    conda:
        ENVS / "bakta.yml"
    threads: 8
    resources:
        mem_mb=double_ram(16 * 1024),
        runtime=2 * 60,
    params:
        work_dir=PROK_ANN / "bakta",
        extra=params["prokaryotes"]["annotate"]["bakta"]["extra"],
        prefix=lambda w: w.mag_id,
    shell:
        """
        bakta \
            --db {input.db} \
            --output {params.work_dir} \
            --prefix {wildcards.mag_id} \
            --threads {threads} \
            --force \
            --keep-contig-headers \
            {params.extra} \
            {input.fasta} \
        2> {log} 1>&2

        for extension in embl faa ffn fna gbff gff3 json tsv; do
            find \
                {params.work_dir} \
                -name "{params.prefix}*.$extension" \
                -exec bgzip {{}} \;
        done 2>> {log} 1>&2
        """


rule prokaryotes__annotate__bakta__all:
    input:
        collect_bakta_files(".tsv.gz"),
