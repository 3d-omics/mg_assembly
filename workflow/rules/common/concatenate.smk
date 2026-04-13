rule concatenate__gzip_text_files:
    """Concatenate multiple text files and output a single text file"""
    input:
        ["file1.fa", "file2.fa.gz"],
    output:
        ["fasta_out.fa.gz"],
    log:
        "fasta_out.log",
    conda:
        ENVS / "concatenate.yml"
    threads: 8
    shell:
        """
        touch {output}

        ( for file in {input} ; do
            if [[ "$file" == *.gz ]] ; then
                gzip \
                    --decompress \
                    --stdout \
                    $file \
            else
                cat $file
            fi \
        done \
        | bgzip \
            --threads {threads} \
        >> {output} \
        ) 2> {log}
        """
