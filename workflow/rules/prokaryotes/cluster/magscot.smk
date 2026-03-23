rule prokaryotes__cluster__magscot__prodigal:
    """Run prodigal over a single assembly"""
    input:
        assembly=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
    output:
        proteins=PROK_MAGSCOT / "{assembly_id}" / "prodigal.faa",
    log:
        PROK_MAGSCOT / "{assembly_id}" / "prodigal.log",
    conda:
        "../../../environments/magscot.yml"
    resources:
        attempt=get_attempt,
    retries: 5
    threads: 24
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=24 * 60,
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.assembly} \
        | parallel \
            --jobs {threads} \
            --block 1M \
            --recstart '>' \
            --pipe \
            --keep-order \
            prodigal \
                -p meta \
                -a /dev/stdout \
                -d /dev/null  \
                -o /dev/null \
        > {output.proteins} \
        ) 2> {log}.{resources.attempt}

        mv {log}.{resources.attempt} {log}
        """


rule prokaryotes__cluster__magscot__hmmsearch_pfam:
    """Run hmmsearch over the predicted proteins of an assembly using Pfam as database

    Note: hmmsearch must be decompressed
    """
    input:
        proteins=PROK_MAGSCOT / "{assembly_id}" / "prodigal.faa",
        hmm=features["magscot"]["pfam_hmm"],
    output:
        PROK_MAGSCOT / "{assembly_id}" / "pfam.tblout.gz",
    log:
        PROK_MAGSCOT / "{assembly_id}" / "pfam.log",
    conda:
        "../../../environments/magscot.yml"
    threads: 4
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=24 * 60,
    shell:
        """
        hmmsearch \
            -o /dev/null \
            --tblout >(pigz --best --processes {threads} > {output}) \
            --noali \
            --notextw \
            --cut_nc \
            --cpu {threads} \
            {input.hmm} \
            {input.proteins} \
        2> {log} 1>&2
        """


rule prokaryotes__cluster__magscot__hmmsearch_tigr:
    """Run hmmsearch over the predicted proteins of an assembly using TIGR as database"""
    input:
        proteins=PROK_MAGSCOT / "{assembly_id}" / "prodigal.faa",
        hmm=features["magscot"]["tigr_hmm"],
    output:
        tblout=PROK_MAGSCOT / "{assembly_id}" / "tigr.tblout.gz",
    log:
        PROK_MAGSCOT / "{assembly_id}" / "tigr.log",
    conda:
        "../../../environments/magscot.yml"
    threads: 4
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=24 * 60,
    shell:
        """
        hmmsearch \
            -o /dev/null \
            --tblout >(pigz --best --processes {threads} > {output.tblout}) \
            --noali \
            --notextw \
            --cut_nc \
            --cpu {threads} \
            {input.hmm} \
            {input.proteins} \
        2> {log} 1>&2
        """


rule prokaryotes__cluster__magscot__join_hmms:
    """Join the results of hmmsearch over TIGR and Pfam

    Note: "|| true" is used to avoid grep returning an error code when no lines are found
    """
    input:
        tigr_tblout=PROK_MAGSCOT / "{assembly_id}" / "tigr.tblout.gz",
        pfam_tblout=PROK_MAGSCOT / "{assembly_id}" / "pfam.tblout.gz",
    output:
        merged=PROK_MAGSCOT / "{assembly_id}" / "hmm.tblout",
    log:
        PROK_MAGSCOT / "{assembly_id}" / "hmm.log",
    conda:
        "../../../environments/magscot.yml"
    shell:
        """
        ( (zgrep -v "^#" {input.tigr_tblout} || true) \
        | awk '{{print $1 "\\t" $3 "\\t" $5}}' ) \
        >  {output.merged} 2>  {log}

        ( (zgrep -v "^#" {input.pfam_tblout} || true) \
        | awk '{{print $1 "\\t" $4 "\\t" $5}}' ) \
        >> {output.merged} 2>> {log}
        """


rule prokaryotes__cluster__magscot__merge_contig_to_bin:
    """Merge the contig to bin files from CONCOCT, MaxBin2 and MetaBAT2

    The output file should have the following format:
    BIN_ID <TAB> CONTIG_ID <TAB> METHOD
    """
    input:
        concoct=PROK_CONCOCT / "{assembly_id}",
        maxbin2=PROK_MAXBIN2 / "{assembly_id}",
        metabat2=PROK_METABAT2 / "{assembly_id}",
    output:
        PROK_MAGSCOT / "{assembly_id}" / "contigs_to_bin.tsv",
    log:
        PROK_MAGSCOT / "{assembly_id}" / "contigs_to_bin.log",
    conda:
        "../../../environments/magscot.yml"
    shell:
        """
        for file in $(find {input.concoct} -name "*.fa.gz" -type f) ; do
            bin_id=$(basename $file .fa)
            zgrep ^">" $file | tr -d ">" \
            | awk -v bin_id=$bin_id '{{print "bin_" bin_id "\\t" $1 "\\tconcoct"}}'
        done > {output} 2> {log}

        for file in $(find {input.maxbin2} -name "*.fa.gz" -type f) ; do
            bin_id=$(basename $file .fa)
            zgrep ^">" $file | tr -d ">" \
            | awk -v bin_id=$bin_id '{{print "bin_" bin_id "\\t" $1 "\\tmaxbin2"}}'
        done >> {output} 2>> {log}

        for file in $(find {input.metabat2} -name "*.fa.gz" -type f) ; do
            bin_id=$(basename $file .fa)
            zgrep ^">" $file | tr -d ">" \
            | awk -v bin_id=$bin_id '{{print "bin_" bin_id "\\t" $1 "\\tmetabat2"}}'
        done >> {output} 2>> {log}
        """


rule prokaryotes__cluster__magscot__run:
    """Run PROK_MAGSCOT over one assembly"""
    input:
        contigs_to_bin=PROK_MAGSCOT / "{assembly_id}" / "contigs_to_bin.tsv",
        hmm=PROK_MAGSCOT / "{assembly_id}" / "hmm.tblout",
    output:
        ar53=PROK_MAGSCOT / "{assembly_id}" / "magscot.gtdb_rel207_ar53.out",
        bac120=PROK_MAGSCOT / "{assembly_id}" / "magscot.gtdb_rel207_bac120.out",
        refined_contig_to_bin=PROK_MAGSCOT
        / "{assembly_id}"
        / "magscot.refined.contig_to_bin.out",
        refined_out=PROK_MAGSCOT / "{assembly_id}" / "magscot.refined.out",
        scores=PROK_MAGSCOT / "{assembly_id}" / "magscot.scores.out",
    log:
        PROK_MAGSCOT / "{assembly_id}/magscot.log",
    conda:
        "../../../environments/magscot.yml"
    params:
        out_prefix=lambda w: PROK_MAGSCOT / w.assembly_id / "magscot",
    resources:
        mem_mb=8 * 1024,
        runtime=12 * 60,
    shell:
        """
        Rscript --vanilla workflow/scripts/MAGScoT/MAGScoT.R \
            --input {input.contigs_to_bin} \
            --hmm {input.hmm} \
            --out {params.out_prefix} \
        2> {log} 1>&2
        """


rule prokaryotes__cluster__magscot__reformat:
    """Reformat the results from PROK_MAGSCOT"""
    input:
        refined_contig_to_bin=PROK_MAGSCOT
        / "{assembly_id}"
        / "magscot.refined.contig_to_bin.out",
    output:
        clean=PROK_MAGSCOT / "{assembly_id}" / "magscot.reformat.tsv",
    log:
        PROK_MAGSCOT / "{assembly_id}" / "magscot.reformat.log",
    conda:
        "../../../environments/magscot.yml"
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        Rscript --vanilla --verbose workflow/scripts/clean_magscot_bin_to_contig.R \
            --input-file {input.refined_contig_to_bin} \
            --output-file {output.clean} \
        2> {log} 1>&2
        """


rule prokaryotes__cluster__magscot__rename:
    """Rename the contigs in the assembly to match the assembly and bin names"""
    input:
        assembly=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
        clean=PROK_MAGSCOT / "{assembly_id}" / "magscot.reformat.tsv",
    output:
        fasta=PROK_MAGSCOT / "{assembly_id}.fa.gz",
    log:
        PROK_MAGSCOT / "{assembly_id}" / "magscot.rename.log",
    conda:
        "../../../environments/magscot.yml"
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        ( python workflow/scripts/reformat_fasta_magscot.py \
            <(gzip -dc {input.assembly}) \
            {input.clean} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.fasta} \
        ) 2> {log}
        """


rule prokaryotes__cluster__magscot__all:
    """Run PROK_MAGSCOT over all assemblies"""
    input:
        [PROK_MAGSCOT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
