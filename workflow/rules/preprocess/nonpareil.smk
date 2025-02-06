rule preprocess__nonpareil__run:
    """Run nonpareil over one sample"""
    input:
        PRE_CLEAN / "{sample_id}.{library_id}_1.fq.gz",
    output:
        redund_val=touch(PRE_NONPAREIL / "{sample_id}.{library_id}.npa"),
        mate_distr=touch(PRE_NONPAREIL / "{sample_id}.{library_id}.npc"),
        log=touch(PRE_NONPAREIL / "{sample_id}.{library_id}.npl"),
        redund_sum=touch(PRE_NONPAREIL / "{sample_id}.{library_id}.npo"),
    log:
        PRE_NONPAREIL / "{sample_id}.{library_id}.log",
    params:
        alg="kmer",
        infer_X=True,
        extra="",
    resources:
        mem_mb=8 * 1024,
        runtime=6 * 60,
    wrapper:
        "v5.2.1/bio/nonpareil/infer"


rule preprocess__nonpareil__plot:
    """Export nonpareil results to json for multiqc"""
    input:
        npo=PRE_NONPAREIL / "{sample_id}.{library_id}.npo",
    output:
        json=PRE_NONPAREIL / "{sample_id}.{library_id}.json",
    log:
        PRE_NONPAREIL / "{sample_id}.{library_id}.json.log",
    wrapper:
        "v5.2.1/bio/nonpareil/plot"


rule preprocess__nonpareil__all:
    """Run nonpareil over all samples and produce JSONs for multiqc"""
    input:
        [
            PRE_NONPAREIL / f"{sample_id}.{library_id}.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
