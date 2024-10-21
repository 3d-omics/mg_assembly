checkpoint prokaryotes__annotate__mags:
    """Separate and decompress all mags from all bins"""
    input:
        assemblies=[MAGSCOT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
    output:
        out_dir=temp(directory(MAGS)),
    log:
        MAGS / ".." / "mags.log",
    conda:
        "__environment__.yml"
    shell:
        """
        mkdir --parents {output.out_dir} 2> {log} 1>&2

        ( gzip \
            --decompress \
            --stdout \
            {input.assemblies} \
        | paste - - \
        | tr -d ">" \
        | tr "@" "\t" \
        | awk \
            '{{print ">" $1 "@" $2 "\\n" $3 > "{output.out_dir}/" $1 ".fa" }}' \
        ) >> {log} 2>&1
        """


localrules:
    prokaryotes__annotate__mags,
