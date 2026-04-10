rule viruses__annotate__dramv__setup:
    input:
        dram_db=features["databases"]["dram"],
    output:
        setup=touch(VIR_DRAMV / "setup.done"),
    log:
        VIR_DRAMV / "setup.log",
    conda:
        ENVS / "dram.yml"
    shell:
        """
        DRAM-setup.py set_database_locations \
            --amg_database_loc          {input.dram_db}/amg_database.*.tsv \
            --dbcan_fam_activities_loc  {input.dram_db}/CAZyDB.*.fam-activities.txt \
            --dbcan_loc                 {input.dram_db}/dbCAN-HMMdb-V*.txt \
            --dbcan_subfam_ec_loc       {input.dram_db}/CAZyDB.*.fam.subfam.ec.txt \
            --description_db_loc        {input.dram_db}/description_db.sqlite \
            --etc_module_database_loc   {input.dram_db}/etc_mdoule_database.*.tsv \
            --function_heatmap_form_loc {input.dram_db}/function_heatmap_form.*.tsv \
            --genome_summary_form_loc   {input.dram_db}/genome_summary_form.*.tsv \
            --kofam_hmm_loc             {input.dram_db}/kofam_profiles.hmm \
            --kofam_ko_list_loc         {input.dram_db}/kofam_ko_list.tsv \
            --module_step_form_loc      {input.dram_db}/module_step_form.*.tsv \
            --peptidase_loc             {input.dram_db}/peptidases.*.mmsdb \
            --pfam_hmm_loc              {input.dram_db}/Pfam-A.hmm.dat.gz \
            --pfam_loc                  {input.dram_db}/pfam.mmspro \
            --viral_loc                 {input.dram_db}/refseq_viral.*.mmsdb \
            --vog_annotations_loc       {input.dram_db}/vog_annotations_latest.tsv.gz \
            --vogdb_loc                 {input.dram_db}/vog_latest_hmms.txt \
        2> {log} 1>&2
        """


rule viruses__annotate__dramv__annotate:
    input:
        fasta=VIR_VIRSORTER2 / "{assembly_id}" / "final-viral-combined-for-dramv.fa",
        tsv=VIR_VIRSORTER2 / "{assembly_id}" / "viral-affi-contigs-for-dramv.tab",
        dram_db=features["databases"]["dram"],
        setup=VIR_DRAMV / "setup.done",
    output:
        annotations=VIR_DRAMV / "annotate" / "{assembly_id}" / "annotations.tsv",
        genes_faa=VIR_DRAMV / "annotate" / "{assembly_id}" / "genes.faa",
        genes_fna=VIR_DRAMV / "annotate" / "{assembly_id}" / "genes.fna",
        genes_gff=VIR_DRAMV / "annotate" / "{assembly_id}" / "genes.gff",
        scaffolds_fna=VIR_DRAMV / "annotate" / "{assembly_id}" / "scaffolds.fna",
        genbank=VIR_DRAMV
        / "annotate"
        / "{assembly_id}"
        / "genbank"
        / "final-viral-combined-for-dramv.gbk",
    log:
        VIR_DRAMV / "annotate" / "{assembly_id}.log",
    conda:
        ENVS / "dram.yml"
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    params:
        workdir=lambda w: VIR_DRAMV / "annotate" / w.assembly_id,
    shell:
        """
        rm \
            --recursive \
            --force \
            --verbose \
            {params.workdir} \
        2> {log} 1>&2

        DRAM-v.py annotate \
            --input_fasta {input.fasta} \
            --output_dir {params.workdir} \
            --skip_trnascan \
            --virsorter_affi_contigs {input.tsv} \
        2>> {log} 1>&2
        """


rule viruses__annotate__dramv__annotate__all:
    input:
        [
            VIR_DRAMV / "annotate" / f"{assembly_id}" / "annotations.tsv"
            for assembly_id in ASSEMBLIES
        ],


use rule csvtk__concat as viruses__annotate__dramv__concatenate_annotations_tsv with:
    input:
        [
            VIR_DRAMV / "annotate" / f"{assembly_id}" / "annotations.tsv"
            for assembly_id in ASSEMBLIES
        ],
    output:
        VIR_DRAMV / "annotations.tsv.gz",
    log:
        VIR_DRAMV / "annotate" / "annotations.log",


use rule concatenate__gzip_text_files as viruses__annotate__dramv__concatenate_genes_fna with:
    input:
        [
            VIR_DRAMV / "annotate" / f"{assembly_id}" / "genes.fna"
            for assembly_id in ASSEMBLIES
        ],
    output:
        VIR_ANN / "dram.genes.fna.gz",
    log:
        VIR_ANN / "dram.genes.fna.log",


use rule concatenate__gzip_text_files as viruses__annotate__dramv__concatenate_genes_faa with:
    input:
        [
            VIR_DRAMV / "annotate" / f"{assembly_id}" / "genes.faa"
            for assembly_id in ASSEMBLIES
        ],
    output:
        VIR_ANN / "dram.genes.faa.gz",
    log:
        VIR_ANN / "dram.genes.faa.log",


use rule concatenate__gzip_text_files as viruses__annotate__dramv__concatenate_scaffolds_fna with:
    input:
        [
            VIR_DRAMV / "annotate" / f"{assembly_id}" / "scaffolds.fna"
            for assembly_id in ASSEMBLIES
        ],
    output:
        VIR_ANN / "dram.scaffolds.fna.gz",
    log:
        VIR_ANN / "dram.scaffolds.fna.log",


use rule concatenate__gzip_text_files as viruses__annotate__dramv__concatenate_genes_gff with:
    input:
        [
            VIR_DRAMV / "annotate" / f"{assembly_id}" / "genes.gff"
            for assembly_id in ASSEMBLIES
        ],
    output:
        VIR_ANN / "dram.genes.gff.gz",
    log:
        VIR_ANN / "dram.genes.gff.log",


rule viruses__annotate__dramv__distill:
    input:
        annotations=VIR_DRAMV / "annotations.tsv.gz",
    output:
        amg_summary=VIR_DRAMV / "amg_summary.tsv.gz",
        vmag_stats=VIR_DRAMV / "vMAG_stats.tsv.gz",
        product=VIR_DRAMV / "product.html",
    log:
        VIR_DRAMV / "distill.log",
    conda:
        ENVS / "dram.yml"
    params:
        outdir=VIR_DRAMV,
        workdir=VIR_DRAMV / "tmp",
    shell:
        """
        rm -rfv {params.workdir} 2> {log}

        DRAM-v.py distill \
            --input_file {input.annotations} \
            --output_dir {params.workdir} \
        2>> {log} 1>&2

        mv \
            {params.workdir}/* \
            {params.outdir}/ \
        2>> {log} 1>&2

        bgzip \
            --threads {threads} \
            {params.outdir}/amg_summary.tsv \
            {params.outdir}/vMAG_stats.tsv \
        2>> {log} 1>&2
        """


rule viruses__annotate__dramv__all:
    input:
        [
            VIR_DRAMV / "annotations.tsv.gz",
            VIR_ANN / "dram.genes.fna.gz",
            VIR_ANN / "dram.genes.faa.gz",
            VIR_ANN / "dram.scaffolds.fna.gz",
            VIR_ANN / "dram.genes.gff.gz",
            VIR_DRAMV / "product.html",
        ],
