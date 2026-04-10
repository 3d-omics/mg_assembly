# reference
# REFERENCE = Path("results/reference/")
# HOSTS = REFERENCE / "hosts"

NULL = [
    Path(os.devnull)
]  # this to make the linter shut up about an absolute path to /dev/null

ENVS = Path("workflow/environments").resolve()

RESULTS = Path("results/")

# Preprocess
PRE = RESULTS / "preprocess"
PRE_READS = PRE / "reads"
PRE_HOSTS = PRE / "hosts"
PRE_FASTP = PRE / "fastp"
PRE_BUILD = PRE / "build"
PRE_BOWTIE2 = PRE / "bowtie2"
PRE_CLEAN = PRE / "clean"
PRE_KRAKEN2 = PRE / "kraken2"
PRE_BRACKEN = PRE / "bracken"
PRE_SINGLEM = PRE / "singlem"
PRE_NONPAREIL = PRE / "nonpareil"

# Assemble
ASSEMBLE = RESULTS / "assemble"
ASMB_MEGAHIT = ASSEMBLE / "megahit"
ASMB_BUILD = ASSEMBLE / "build"
ASMB_BOWTIE2 = ASSEMBLE / "bowtie2"
ASMB_COVERM = ASSEMBLE / "coverm"
ASMB_QUAST = ASSEMBLE / "quast"
ASMB_KRAKEN2 = ASSEMBLE / "kraken2"

# Prokaryotes
PROK = RESULTS / "prokaryotes"

## Prokaryotes - Cluster
PROK_CLUSTER = PROK / "cluster"
PROK_CONCOCT = PROK_CLUSTER / "concoct"
PROK_METABAT2 = PROK_CLUSTER / "metabat2"
PROK_MAXBIN2 = PROK_CLUSTER / "maxbin2"
PROK_MAGSCOT = PROK_CLUSTER / "magscot"
PROK_PRODIGAL = PROK_CLUSTER / "prodigal"

## Prokaryotes - Annotate
PROK_ANN = PROK / "annotate"
PROK_MAGS = PROK_ANN / "mags"
PROK_GTDBTK = PROK_ANN / "gtdbtk"
PROK_QUAST = PROK_ANN / "quast"
PROK_DRAM = PROK_ANN / "dram"
PROK_DREP = PROK_CLUSTER / "drep"


## Prokaryotes - Quantify
PROK_QUANT = PROK / "quantify"
PROK_BUILD = PROK_QUANT / "build"
PROK_BOWTIE2 = PROK_QUANT / "bowtie2"
PROK_COVERM = PROK_QUANT / "coverm/"


# Viruses
VIR = RESULTS / "viruses"

## Viruses - Cluster
VIR_CLUSTER = VIR / "cluster"
VIR_GENOMADC = VIR_CLUSTER / "genomad"
# VIR_CHECKVC = VIR_CLUSTER / "checkv"
VIR_DEDUPE = VIR_CLUSTER / "dedupe"
VIR_MMSEQS = VIR_CLUSTER / "mmseqs"

## Viruses - Annotation
VIR_ANN = VIR / "annotate"
VIR_VIRSORTER2 = VIR_ANN / "virsorter2"
VIR_GENOMADA = VIR_ANN / "genomad"
VIR_DRAMV = VIR_ANN / "dramv"
VIR_QUAST = VIR_ANN / "quast"
VIR_CHECKV = VIR_ANN / "checkv"

## Viruses - Quantify
VIR_QUANT = VIR / "quantify"
VIR_BUILD = VIR_QUANT / "build"
VIR_BOWTIE2 = VIR_QUANT / "bowtie2"
VIR_COVERM = VIR_QUANT / "coverm"
