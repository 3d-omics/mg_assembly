def get_attempt(wildcards, attempt):
    """Get the number of attempt in resources"""
    return attempt


def double_ram(initial_memory):
    """Double the memory for each attempt"""
    return lambda wildcards, attempt: initial_memory * 2 ** (attempt - 1)


def collect_mag_outputs(path_pattern):
    """Return a checkpoint input function that collects per-MAG paths.

    path_pattern: Path or str with {mag_id} placeholder(s).
    Example: PROK_ANN / "bakta" / "annotate" / "{mag_id}" / "{mag_id}.faa"
    """

    def _inner(wildcards):
        checkpoints.prokaryotes__annotate__mags.get()
        mag_ids = glob_wildcards(PROK_MAGS / "{mag_id}.fa").mag_id
        return expand(path_pattern, mag_id=mag_ids)

    return _inner
