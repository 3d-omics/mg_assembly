def collect_bakta_files(ext):
    """Collect specific Bakta output files for all MAGs."""
    return collect_mag_outputs(PROK_ANN / "bakta" / f"{{mag_id}}{ext}")
