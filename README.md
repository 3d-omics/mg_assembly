# Snakemake workflow: `Bioinfo_Macro_Genome_Resolved_Metagenomics`

A Snakemake workflow for Genome Resolved Metagenomics

##Set up required softwares

## Usage
  ```
  #Clone the git repository in your terminal
  git clone git@github.com:3d-omics/Bioinfo_Macro_Genome_Resolved_Metagenomics.git
  #Change directory to the one you cloned in the previous step
  cd Bioinfo_Macro_Genome_Resolved_Metagenomics
  #Activate conda environment where you have snakemake
  conda activte Snakemake
  #run the pipeline with the test data, it will download all the necesary software through conda. It should take less than 5 minutes.
  snakemake --use-conda --jobs 8 all
  ```

- Run it with your own data:
  - Edit `config/samples.tsv` and add your samples and where are they located. Here is an example of the tsv table filled with the information


  - Run the pipeline
     ```
     snakemake --use-conda --jobs 8 all
     #(slurm users), there is a script called run_slurm in the cloned directory that you can directly use to launch the pipeline on a slurm         cluster, you can modify the parameters or direclty execute it as it is
     ./run_slurm
     ```

## Features
- FASTQ processing with `fastp`
- mapping of preprocessed reads against the host with `bowtie2`
- assembly of non-host reads with `megahit`
- binning with `metabin`
- coverage evaluation with `coverm`

## DAG

![image](https://github.com/3d-omics/Bioinfo_Macro_Genome_Resolved_Metagenomics/assets/103645443/08147c8a-8939-437f-aa12-2edede627ed6)


## References


