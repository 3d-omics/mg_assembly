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
            if [[ $file =~ \.gz$ ]] ; then
                command="gzip --decompress --stdout"
            else
                command="cat"
            fi

            $command $file
        done \
        | bgzip \
            --threads {threads} \
            --stdout \
        >> {output} \
        ) 2> {log}
        """
