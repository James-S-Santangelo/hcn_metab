## Analysis of RNA-seq data exploring tissue-specific expression of HCN in white clover

### Description of repository

This repository contains code necessary to reproduce the RNA-seq analyses from [this paper]().
Raw read files have been deposited on NCBI (BioProject [PRJNA1124556](https://www.ncbi.nlm.nih.gov/sra/PRJNA1124556). The pipeline uses `Conda`, `Snakemake`, 
and `Singularity` for workflow management and reproducibility. All Snakefiles and directories are well-documented, 
but here is a brief overview of the pipeline and directories in this repository:

1. Align Raw RNA-seq reads to the white clover reference genome (GenBank [GCA_030408175.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030408175.1/)
2. Generate normalized feature counts around HCN loci and perform differential gene expression analysis using DESeq2 

#### Overview of directories

- [config](./config): Snakemake configuration files for different clusters.
- [resources](./resources): Text files used in pipeline (e.g., sample information)
- [workflow](./workflow): Main Snakemake workflow with rules, environments, scripts, etc

### Using the pipeline

This pipeline requires `Conda`. Note the reference genome and raw reads need to be downloaded first,
and their paths should be specified in the config file (see below)

- A minimal installation of `Conda` (i.e., Miniconda) can be installed by following the instructions for your platform [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Assuming `Conda` is installed, the this repository's `Conda` environment can be replicated by running the following command:

`conda env create -f environment.yaml -n hcn_metab`

This will create a `Conda` environment named _hcn\_metab_ containing a minimal set of dependencies required to run the pipeline. This environment additinally contains dependencies to run Jupter Notebooks.

After activating the environment (`conda activate hcn_metab`), the pipeline can be executed from the [workflow](./workflow) directory by running a command that looks something like:

`snakemake --use-conda --configfile ../config/<configfile> --notemp -j <cores>`

for local execution. Here, `<path>` is the path on the cluster from which files will be read/written (e.g., `/scratch`), 
`<configfile>` is one of the configfiles in the [config](./config) directory that needs to be modified to match the paths on your system, 
and `<cores>` is the number of cores available for executing parallel processes. 

