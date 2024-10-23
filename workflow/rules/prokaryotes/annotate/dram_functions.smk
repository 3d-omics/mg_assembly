def collect_dram_annotate(wildcards):
    checkpoint_output = checkpoints.prokaryotes__annotate__mags.get().output[0]
    mag_ids = glob_wildcards(MAGS / "{mag_id}.fa").mag_id
    return [PROK_ANN / "dram.annotate" / mag_id for mag_id in mag_ids]
