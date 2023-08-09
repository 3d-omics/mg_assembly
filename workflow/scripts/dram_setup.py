#!/usr/bin/env python3


import subprocess
import os
import sys


def get_kofam_hmm_loc(folder):
    return folder + "/kofam_profiles.hmm"


def get_kofam_ko_list_loc(folder):
    return folder + "/kofam_ko_list.tsv"


def get_uniref_loc(folder):
    files = os.listdir(folder)
    uniref_file = [
        file for file in files if ("uniref90" in file and file.endswith(".mmsdb"))
    ][0]
    return f"{folder}/{uniref_file}"


def get_pfam_loc(folder):
    return folder + "/pfam.mmspro"


def get_pfam_hmm_loc(folder):
    return folder + "/Pfam-A.hmm.dat.gz"


def get_dbcan_loc(folder):
    return folder + "/dbCAN-HMMdb-V9.txt"


def get_dbcan_fam_activities_loc(folder):
    files = os.listdir(folder)
    dbcan_fam_activities_file = [
        file
        for file in files
        if ("CAZyDB" in file and file.endswith("fam-activities.txt"))
    ][0]
    return f"{folder}/{dbcan_fam_activities_file}"


def get_viral_loc(folder):
    files = os.listdir(folder)
    viral_file = [
        file for file in files if ("refseq_viral" in file and file.endswith(".mmsdb"))
    ][0]
    return f"{folder}/{viral_file}"


def get_peptidase_loc(folder):
    files = os.listdir(folder)
    peptidase_file = [
        file for file in files if ("peptidases" in file and file.endswith(".mmsdb"))
    ][0]
    return f"{folder}/{peptidase_file}"


def get_vogdb_loc(folder):
    return folder + "/vog_latest_hmms.txt"


def get_vog_annotations_loc(folder):
    return folder + "/vog_annotations_latest.tsv.gz"


def get_genome_summary_form_loc(folder):
    files = os.listdir(folder)
    genome_summary_form_file = [
        file
        for file in files
        if ("genome_summary_form" in file and file.endswith(".tsv"))
    ][0]
    return f"{folder}/{genome_summary_form_file}"


def get_module_step_form_loc(folder):
    files = os.listdir(folder)
    module_step_form_file = [
        file for file in files if ("module_step_form" in file and file.endswith(".tsv"))
    ][0]
    return f"{folder}/{module_step_form_file}"


def get_etc_module_database_loc(folder):
    files = os.listdir(folder)
    etc_module_database_file = [
        file
        for file in files
        if ("etc_mdoule_database" in file and file.endswith(".tsv"))
    ][0]
    return f"{folder}/{etc_module_database_file}"


def get_function_heatmap_form_loc(folder):
    files = os.listdir(folder)
    function_heatmap_form_file = [
        file
        for file in files
        if ("function_heatmap_form" in file and file.endswith(".tsv"))
    ][0]
    return f"{folder}/{function_heatmap_form_file}"


def get_amg_database_loc(folder):
    files = os.listdir(folder)
    amg_database_file = [
        file for file in files if ("amg_database" in file and file.endswith(".tsv"))
    ][0]
    return f"{folder}/{amg_database_file}"


def get_description_db_loc(folder):
    return folder + "/description_db.sqlite"


def compose_dram_setup_set_database_locations(folder):
    command = [
        "DRAM-setup.py",
        "set_database_locations",
        "--kofam_hmm_loc",
        get_kofam_hmm_loc(folder),
        "--kofam_ko_list_loc",
        get_kofam_ko_list_loc(folder),
        "--uniref_loc",
        get_uniref_loc(folder),
        "--pfam_loc",
        get_pfam_loc(folder),
        "--pfam_hmm_loc",
        get_pfam_hmm_loc(folder),
        "--dbcan_loc",
        get_dbcan_loc(folder),
        "--dbcan_fam_activities_loc",
        get_dbcan_fam_activities_loc(folder),
        "--vogdb_loc",
        get_vogdb_loc(folder),
        "--vog_annotations_loc",
        get_vog_annotations_loc(folder),
        "--viral_loc",
        get_viral_loc(folder),
        "--peptidase_loc",
        get_peptidase_loc(folder),
        "--description_db_loc",
        get_description_db_loc(folder),
        "--genome_summary_form_loc",
        get_genome_summary_form_loc(folder),
        "--module_step_form_loc",
        get_module_step_form_loc(folder),
        "--etc_module_database_loc",
        get_etc_module_database_loc(folder),
        "--function_heatmap_form_loc",
        get_function_heatmap_form_loc(folder),
        "--amg_database_loc",
        get_amg_database_loc(folder),
    ]
    return " ".join(command)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 dram_setup.py <folder>")
        exit(1)
    folder = sys.argv[1]
    command = compose_dram_setup_set_database_locations(folder)
    subprocess.run(command, shell=True)
