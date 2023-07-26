# EpicTope

Simplified code for identifying ideal epitope tag insertion sites for proteins.

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)

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
## Installation

Once dependencies have been installed, you will need to install the `epictope` R package from this github repo. We can use the `remotes` package to do this.

```
install.packages("remotes")
remotes::install_github("henrichung/epitope_tag")
```
## Usage


#### Workflows
Example workflows with the `epictope` package are available in the **vignettes** folder. Workflows are available as both R Markdown Documents and Jupyter notebooks. These workflows describe the `epictope` workflow step by step in an interactive session or an IDE.

#### Wrappers
Alternatively, the wrapper scripts `install.R` and `single_score.R` are provided in the **scripts** folder of this repo to enable one-touch instant operation.
To run, download the `install.R` and `single_score.R` scripts and place them into your current project directory in a *code* folder. Your project directory should look as follows;

Folder Structure
- `project/`: Project folder
  - `code/`: 
    - `install.R`: downloads and creates blasteable databases
    - `single_score.R`: runs tagging software on a single protein


From the terminal, the scripts can be run as follows.
```
Rscript code/install.R
Rscript code/single_score "P57102" # replace 'P57102' with the UniprotID for your protein of interest.
```

Additionally, an example `config.R` file is provided to change any of the user-customizeable values during feature scoring. The values listed in the example config are the default values used by `epictope`.
