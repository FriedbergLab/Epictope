# EpicTope

Simplified code for identifying ideal epitope tag insertion sites for proteins.

## Table of Contents

- [Config](#config)
- [Installation](#installation)
- [Dependencies](#dependencies)

## Config

A user customizeable config file determines the folder structure, species for multiple sequence alignment, feature weights, and query input protein. 

## Installation

A installation script installs necessary R packages, downloads CDS files from NCBI FTP, and creates protein blastable databases.

```
# Rscript install.R
```
## Dependencies

To calculate the multiple sequence alignment and secondary characteristics, `epictope` relies on local installs of BLAST, muscle, and dssp. Fortunately, these packages can be installed using conda, an open-source package management system and environment management system that runs on Windows, macOS, and Linux. Conda installers can be found at the Anaconda[website](https://www.anaconda.com/). Once installed, you may run the follow commands to install the requisite packages. These commands will create a conda environment named "epictope", and install the requisite packages into that environment. 

```
conda create -n epictope
conda activate epictope
conda install -c bioconda blast
conda install -c bioconda muscle
conda install -c salilab dssp
conda install -c conda-forge r-base
```
#### Installation

Once installed, you will need to open R and install the `remotes` packages. This will allow you to install the `epictope` package from this Github repo.

```
install.packages("remotes")
remotes::install_github("henrichung/epitope_tag")
```


## Usage

Example workflows with the `epictope` package are available in the "vignettes" folder. Workflows are available as both R Markdown Documents and Jupyter notebooks. These 

```
# Rscript single_score.R
```

Folder Structure
- `R/`: R implementation
- `python`: Python implementation 
- `data/`: 
  - `models/`: Folder for downloaded query pdb files
  - `CDS`: Folder for CDS files for species in MSA
- `outputs/`: 
