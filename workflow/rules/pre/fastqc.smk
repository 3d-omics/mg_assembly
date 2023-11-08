rule pre_fastp_fastqc:
    """Run fastqc over all libraries after fastp"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in "1 2".split(" ")
            for extension in "html zip".split(" ")
        ],


rule pre_nonhost_fastqc:
    """Run fastqc over all libraries after fastp"""
    input:
        [
            PRE_BOWTIE2
            / f"non{genome}/{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in "1 2".split(" ")
            for extension in "html zip".split(" ")
            for genome in HOST_NAMES
        ],
