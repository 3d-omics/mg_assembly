rule magscot_prodigal_one:
    input:
        assembly=ASSEMBLY_RENAME / "{assembly_id}.fa",
    output:
        proteins=MAGSCOT / "{assembly_id}/prodigal.faa",
        genes=MAGSCOT / "{assembly_id}/prodigal.ffn",
    log:
        MAGSCOT / "{assembly_id}/prodigal.log",
    conda:
        "../../envs/metabin.yml"
    threads: 24
    shell:
        """
        (cat {input.assembly} \
        | parallel \
            --jobs {threads} \
            --block 1M \
            --recstart '>' \
            --pipe \
            prodigal \
                -p meta \
                -a {output.proteins}.{{#}}.faa \
                -d {output.genes}.{{#}}.ffn \
                -o /dev/null \
        ) 2> {log} 1>&2
        cat {output.proteins}.*.faa > {output.proteins} 2>> {log}
        cat {output.genes}.*.ffn > {output.genes} 2>> {log}
        rm -f {output.proteins}.*.faa {output.genes}.*.ffn 2>> {log} 2>&1
        """


rule magscot_hmmsearch_pfam_one:
    input:
        proteins=MAGSCOT / "{assembly_id}/prodigal.faa",
        hmm=features["magscot"]["pfam_hmm"],
    output:
        tblout=MAGSCOT / "{assembly_id}/pfam.tblout",
    log:
        MAGSCOT / "{assembly_id}/pfam.log",
    conda:
        "../../envs/metabin.yml"
    threads: 2
    shell:
        """
        hmmsearch \
            -o /dev/null \
            --tblout {output.tblout} \
            --noali \
            --notextw \
            --cut_nc \
            --cpu {threads} \
            resources/magscot/gtdbtk_rel207_tigrfam.hmm \
            {input.proteins} \
        2> {log} 1>&2
        """


rule magscot_hmmsearch_tigr_one:
    input:
        proteins=MAGSCOT / "{assembly_id}/prodigal.faa",
        hmm=features["magscot"]["tigr_hmm"],
    output:
        tblout=MAGSCOT / "{assembly_id}/tigr.tblout",
    log:
        MAGSCOT / "{assembly_id}/tigr.log",
    conda:
        "../../envs/metabin.yml"
    threads: 2
    shell:
        """
        hmmsearch \
            -o /dev/null \
            --tblout {output.tblout} \
            --noali \
            --notextw \
            --cut_nc \
            --cpu {threads} \
            resources/magscot/gtdbtk_rel207_tigrfam.hmm \
            {input.proteins} \
        2> {log} 1>&2
        """


rule magscot_join_hmm_one:
    input:
        tigr_out=MAGSCOT / "{assembly_id}/tigr.tblout",
        pfam_out=MAGSCOT / "{assembly_id}/pfam.tblout",
    output:
        merged=MAGSCOT / "{assembly_id}/hmm.tblout",
    log:
        MAGSCOT / "{assembly_id}/hmm.log",
    conda:
        "../../envs/metabin.yml"
    shell:
        """
        grep -v "^#" {input.tigr_out} | awk '{{print $1"\\t"$3"\\t"$5}}' >  {output.merged} 2>  {log}
        grep -v "^#" {input.pfam_out} | awk '{{print $1"\\t"$4"\\t"$5}}' >> {output.merged} 2>> {log}
        """


rule magscot_compose_contig_to_bin_concoct_one:
    input:
        BIN_METAWRAP / "{assembly_id}/concoct_bins",
    output:
        MAGSCOT / "{assembly_id}/concoct.contigs_to_bin.tsv",
    log:
        MAGSCOT / "{assembly_id}/concoct.contigs_to_bin.log",
    conda:
        "../../envs/metabin.yml"
    shell:
        """
        (grep -H ^">" {input}/*.fa \
        | parallel -j 1 echo {{/}} \
        | sed 's/\.fa:>/\\t/' \
        | awk '{{print $0"\\tconcoct"}}' \
        > {output} \
        ) 2> {log}
        """


rule magscot_compose_contig_to_bin_maxbin2_one:
    input:
        BIN_METAWRAP / "{assembly_id}/maxbin2_bins",
    output:
        MAGSCOT / "{assembly_id}/maxbin2.contigs_to_bin.tsv",
    log:
        MAGSCOT / "{assembly_id}/maxbin2.contigs_to_bin.log",
    conda:
        "../../envs/metabin.yml"
    shell:
        """
        (grep -H ^">" {input}/*.fa \
        | parallel -j 1 echo {{/}} \
        | sed 's/\.fa:>/\\t/' \
        | awk '{{print $0"\\tmaxbin2"}}' \
        > {output} \
        ) 2> {log}
        """


rule magscot_compose_contig_to_bin_metabat2_one:
    input:
        BIN_METAWRAP / "{assembly_id}/metabat2_bins",
    output:
        MAGSCOT / "{assembly_id}/metabat2.contigs_to_bin.tsv",
    log:
        MAGSCOT / "{assembly_id}/metabat2.contigs_to_bin.log",
    conda:
        "../../envs/metabin.yml"
    shell:
        """
        (grep -H ^">" {input}/*.fa \
        | parallel -j 1 echo {{/}} \
        | sed 's/\.fa:>/\\t/' \
        | awk '{{print $0"\\tmetabat2"}}' \
        > {output} \
        ) 2> {log}
        """


rule magscot_merge_contig_to_bin_one:
    input:
        MAGSCOT / "{assembly_id}/concoct.contigs_to_bin.tsv",
        MAGSCOT / "{assembly_id}/maxbin2.contigs_to_bin.tsv",
        MAGSCOT / "{assembly_id}/metabat2.contigs_to_bin.tsv",
    output:
        MAGSCOT / "{assembly_id}/contigs_to_bin.tsv",
    log:
        MAGSCOT / "{assembly_id}/contigs_to_bin.log",
    conda:
        "../../envs/metabin.yml"
    shell:
        """
        cat {input} > {output} 2> {log}
        """


rule magscot_run_one:
    input:
        contigs_to_bin=MAGSCOT / "{assembly_id}/contigs_to_bin.tsv",
        hmm=MAGSCOT / "{assembly_id}/hmm.tblout",
    output:
        ar53=MAGSCOT / "{assembly_id}/magscot.gtdb_rel207_ar53.out",
        bac120=MAGSCOT / "{assembly_id}/magscot.gtdb_rel207_bac120.out",
        refined_contig_to_bin=MAGSCOT
        / "{assembly_id}/magscot.refined.contig_to_bin.out",
        refined_out=MAGSCOT / "{assembly_id}/magscot.refined.out",
        scores=MAGSCOT / "{assembly_id}/magscot.scores.out",
    log:
        MAGSCOT / "{assembly_id}/magscot.log",
    conda:
        "../../envs/metabin.yml"
    params:
        out_prefix=lambda wildcards: MAGSCOT / f"{wildcards.assembly_id}/magscot",
    shell:
        """
        Rscript --vanilla workflow/scripts/MAGScoT/MAGScoT.R \
            --input {input.contigs_to_bin} \
            --hmm {input.hmm} \
            --out {params.out_prefix} \
        2> {log} 1>&2
        """


rule magscot_reformat_one:
    input:
        refined_contig_to_bin=MAGSCOT
        / "{assembly_id}/magscot.refined.contig_to_bin.out",
    output:
        clean=MAGSCOT / "{assembly_id}/magscot.reformat.tsv",
    log:
        MAGSCOT / "{assembly_id}/magscot.reformat.log",
    conda:
        "../../envs/metabin.yml"
    shell:
        """
        Rscript --vanilla workflow/scripts/clean_magscot_bin_to_contig.R \
            --input-file {input.refined_contig_to_bin} \
            --output-file {output.clean} \
        2> {log} 1>&2
        """


rule magscot_rename_one:
    input:
        assembly=ASSEMBLY_RENAME / "{assembly_id}.fa",
        clean=MAGSCOT / "{assembly_id}/magscot.reformat.tsv",
    output:
        fasta=MAGSCOT / "{assembly_id}.fa",
    log:
        MAGSCOT / "{assembly_id}/magscot.refined.log",
    conda:
        "../../envs/metabin.yml"
    shell:
        """
        python workflow/scripts/reformat_fasta_magscot.py \
            {input.assembly} \
            {input.clean} \
        > {output.fasta} 2> {log}
        """


rule magscot_split_into_bins:
    input:
        fasta=MAGSCOT / "{assembly_id}.fa",
    output:
        bins=directory(MAGSCOT / "{assembly_id}/bins"),
    log:
        MAGSCOT / "{assembly_id}/bins.log",
    conda:
        "../../envs/metabin.yml"
    shell:
        """
        mkdir -p {output.bins} 2> {log}
        (seqtk seq results/binning/magscot/all.fa \
        | paste - -  \
        | tr "." "\\t" \
        | awk '{{print $1"."$2"."$3"\\n"$4 > "{output.bins}/"$2".fa"}}'
        ) 2> {log}
        """


rule magscot:
    input:
        [MAGSCOT / f"{assembly_id}/bins" for assembly_id in ASSEMBLIES],
