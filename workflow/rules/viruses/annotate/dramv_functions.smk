def collect_dramv_annotate(wildcards):
    checkpoint_output = checkpoints.viruses__annotate__dramv__contigs.get().output[0]
    contig_ids = glob_wildcards(VIR_DRAMV / "contigs" / "{contig_id}.fa").contig_id
    return [VIR_DRAMV / "annotate" / contig_id for contig_id in contig_ids]
