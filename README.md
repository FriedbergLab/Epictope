# Epictope

Software for predicting epitope tag insertion sites in proteins.

Epictope is an R pipeline to identify epitope tag insertion sites for proteins of interest. It uses four features of protein structure; sequence conservation, secondary structure, disordered binding regions, and relative solvent accessabilty to predict suitable internal locations for tag insertion.

This repository contains the code source of the R Epictope package, step-by-step R Markdown and Jupyter notebooks to run the complete workflow, wrapped scripts for simplified workflow executation, and instructions to adjust the weight and effect of each considered feature. The package requires local installations of BLAST, MUSCLE, and DSSP to run (those can b eisntalled as a apckage with EpicTope). You will need at elast 3GB of disk space.

## Table of Contents
- [Methodology](#methodology)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)

## Methodology

### Sequence conservation

Sequence conservation is used to guide internal epitope-tagging approaches. Regions of relatively low conservation are unlikely to be involved in the critical function of the protein. To identify these regions for a protein of interest, we first BLAST the query protein against the proteomes of a diverse set of model organisms. By default, we compare the query sequence against the proteomes of _Mus musculus_ (mouse), _Bos taurus_ (cow), _Canis lupus familiaris_ (dog), Gallus gallus (chicken), _Homo sapiens_ (human), _Takifugu rubripes_ (pufferfish), and _Xenopus tropicalis_ (western clawed frog). Using BLAST, we identify the highest scoring match in each organism, sorted by the lowest E-value. We then align the retrieved sequences with the query protein using MUSCLE, a multiple sequence alignment program, and calculate the shannon entropy at each position. We use [Shannon entropy](https://en.wikipedia.org/wiki/Entropy_(information_theory)) as a simple measure of the calculate the variability of amino acids at each position in the alignment. A lower Shannon entropy indicates low variability, or high sequence conservation at the position, and it should therefore be avoided for tag insertion. Conversly, a high Shannon entropy indicates a relatively low degree of sequence conservation, and potential suitability for tagging.


<figure style="display: inline-block; text-align: center;">
  <img src="images/msa.png" alt="Alt text" title="Tcf21 Multiple Sequence Alignment." width="75%">
  <figcaption>Example Multiple Sequence Alignment for Tcf21 protein sequences from position 1 to 70. Amino acids identical between all species are in red, identical between at least three out of five species in blue. From this alignment, red regions would be extremely unfavorable to tag insertion. Figure from https://doi.org/10.1038/srep36986</figcaption>. Reproduced under Creative Commons CC-BY.
</figure>

### Solvent accessibility

 Relative Solvent Accessibility (RSA) is a measure of the surface area of a folded protein that is accessible to a solvent, typically the cytoplasmic fluid. It is calculated by dividing the solvent accessible surface area (SASA) of an amino acid by the maximum possible solvent accessible surface area for that residue. SASA values are assigned with Define Secondary Structure of Proteins (DSSP). The DSSP program defines secondary structure, geometrical features and solvent exposure of proteins, given atomic coordinates in Protein Data Bank (PDB) format. Values used for the maximum possible solvent accessible surface area were taken from [this study](https://doi.org/10.1371/journal.pone.0080635). We use the [Alphafold2 predicted structure from the European Bioinforamtics Institute (EBI)](https://alphafold.ebi.ac.uk/) as the source PDB for DSSP calculations.

Insert _non copyrighted_ figure here.



<!--
Insert non copyrighted figure here.
<figure style="display: inline-block; text-align: center;">
  <img src="images/rsa.png" alt="Alt text" title="Solvent Accessible Regions." width="75%">
  <figcaption> Crystal structures of proteins with varying relative solvent-accessible surface area (given in parentheses). Buried residues, or inaccessible residues, are colored blue and solvent exposed residues are colored red. Figure from https://doi.org/10.1016/j.jmb.2013.06.019 .</figcaption>
</figure>
-->

### Secondary structure

Secondary structure is the local spatial conformation of the polypeptide backbone for the protein of interest. Certain structures, such as alpha helices or beta sheets, are more defined and disruption of these structure is likely to affect protein structure. As with solvent accessibility, we use DSSP to define the secondary structure of the protein from its PDB file. By default, we assign helices (GHI) and sheets (E) feature scores of 0. Hydrogen bonded turns (T), residues in isolated Beta bridges (B), and bends (S) scores of 0.5, and coils scores of 1. For all features, higher values indicate greater suitability for tag insertion. 

### Disordered binding 

Disordered binding regions are sections of a protein that do not have a well-defined structure on their own, but can undergo a disorder-to-order transition when they bind to specific protein partners. To avoid these regions, we use [ANCHOR2](https://iupred2a.elte.hu/), a tool that analyzes an amino acid sequence and returns a score of intrinsic disorder depending on a model of the estimated energy potential for residue interactions. To maintain consistency with other features, the disordered binding feature score is taken as 1 minus the ANCHOR2 score.


## Installation

### Dependencies
To calculate the multiple sequence alignment and secondary characteristics, Epictope relies on local installs of BLAST, muscle, and DSSP. These packages can be installed using conda, an open-source package management system and environment management system that runs on Windows, macOS, and Linux. Conda installers can be found at the Anaconda [website](https://www.anaconda.com/). Once installed, run the following commands to install the requisite packages. These commands will create a conda environment named "epictope", and install the requisite packages into that environment.

=======
## Requirements

Installing `epictope` and its dependencies will require at least 3Gb of disk space. Users should also be familar with using conda, a package manager for macOS/linux and Windows. Conda does not need to be used if users already have access to  installations of BLAST, muscle, and dssp, either locally or on an HPC environment. For users familiar with R, `epictope` can be run interactively through an R session or with an IDE such as RStudio. 

## Dependencies


To calculate the multiple sequence alignment and secondary characteristics, `epictope` relies on local installs of BLAST, muscle, and dssp. These packages can be installed using conda, an open-source package management system and environment management system that runs on Windows, macOS, and Linux. Conda installers can be found at the Anaconda [website](https://www.anaconda.com/). Once installed, you may run the follow commands to install the requisite packages. These commands will create a conda environment named "epictope", and install the requisite packages into that environment. 

For macOS/Linux, commands are issued at the terminal. For Windows, it is recommended to use conda through the Anaconda Prompt. Detailed instructions for Windows can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)

Dependencies can be installed using the provided `epictope_environment.yml` using the following commands.
```
conda env create --file=epictope_environment.yml
conda activate epictope
```

Alternatively, dependencies can be installed separately by first creating an empty environment, activating the environment, and issuing install commands individually. This can resolve conflict issues that sometimes occur when installing via yml.

```
conda create -n epictope
conda activate epictope
conda install -c bioconda blast
conda install -c bioconda muscle
conda install -c salilab dssp
conda install -c conda-forge r-base
conda install -c conda-forge r-stringi
conda install -c conda-forge r-openssl
```

We also provide a simple wrapper scripts `environment_install.sh` for macOS/Linux and `environment_install.bat` for Windows, which issues the individual conda install commands as a single script. 

```
# for macOS/Linux
chmod +x environment_install.sh
./environment_install.sh

# for Windows
C:\ProgramData\Anaconda3\Scripts\Activate
environment_install.bat
```

#### Example dependency installations on Linux/macOS

Method 1
```
conda env create --file=epictope_environment.yml
conda activate epictope
```

Method 2
```
chmod +x environment_install.sh
./environment_install.sh
```

#### Example dependency installations on Windows

Method 1
```
C:\ProgramData\Anaconda3\Scripts\Activate  
conda env create --file=epictope_environment.yml
conda activate epictope
```

Method 2
```
C:\ProgramData\Anaconda3\Scripts\Activate
environment_install.bat
```
### Installing EpicTope 

Epictope is distributed as a R package. You can install it from this github repository using the *install_github* function from eiher the `remotes` or `devtools` R packages. This function is equivalent to the *install.packages()* function in base R, and should be entered in R through an interactive session or an IDE like RStudio.

```
# using remotes
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("henrichung/epitope_tag")

# using devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("henrichung/epitope_tag")
```

We also provide additional resources outside the main packages, such as scripts to automate conda installation, wrapper scripts for the epictope workflow, and detailed step-by-step workflows as R Markdown Documents and Jupyter Notebooks. These can be downloaded from this repository through the webpage, or using git. Git is available on both MacOS/Linux and Windows machines

```
git clone https://github.com/henrichung/epitope_tag
```

## Usage

#### Workflow notebooks
Example workflows with the `epictope` package are available in the **vignettes** folder. Workflows are available as both [R Markdown Documents](https://rmarkdown.rstudio.com/) and [Jupyter](https://jupyter.org/) notebooks. These workflows go through the `epictope` workflow step by step in an interactive session or an IDE. 

#### Macro scripts
Alternatively, the scripts `install.R` and `single_score.R` are provided in the **scripts** folder of this repo to enable one-command operation.
To run, download the `install.R` and `single_score.R` scripts from this repository either directly from the github page or using git clone.

- [install.R](https://github.com/henrichung/epitope_tag/blob/main/scripts/install.R) - This script downloads the proteomes for the species used in the multiple sequence alignment. It then converts these sequences in useable files for BLAST. The script checks for and, if not existent, creates data folders in the current working directory to store these files. These file needs to be re-run if the user changes the species considered in the multiple sequence alignment. 
- [single_score.R](https://github.com/henrichung/epitope_tag/blob/main/scripts/single_score.R) - This script takes a UniprotID as input and performs the epictope workflow for that protein. It retrieves the amino acid sequence and Alphafold2 predicted structure for the protein. It then BLASTs the protein against the proteomes of the animals used in the multiple sequence alignment, retrieves the highest scoring match (score measued by lowest E-value), and aligns the matched proteins along with the query in a multiple sequence alignment. It separately determines the secondary structure, solvent accessibility, and disordered binding regions for the protein, and finally combines all feature scores into a dataframe. The dataframe annotates each residue position with its feature scores and final tagging score. This file is saved to an /outputs folder with the name of the protein followed by '_score.csv'. For example, the protein used in the example below save a "outputs/P57102_score.csv" file.

From the terminal, these scripts can be run as follows.
```
Rscript install.R 
Rscript single_score.R "P57102" # replace 'P57102' with the UniprotID for your protein of interest.
```

Each script can also be opened in an IDE such as Rstudio, and run interactively line by line.

#### User configuration

The scoring function used by `epictope` sums the calculated scores for each feature, assigning equal weight to thesecondary structure, disordered binding regions, and solvent accessabilty. Sequence conservation is weighted slightly higher, at 1.5 times the other features. The weight for each feature can be adjusted by the user using a "config.R" file. This file (shorthand for configuration) is used to adjust some of the tuneable parameters in epictope. The configuration file allows the user to define the species used in the multiple sequence alignment, the values used to score the tag suitability of secondary structures, and the maximum solvent accessibility values to determine solvent accessibility. By default, `epictope` will look for a config.R file in the working directory. If a file is not found, it will use default values. The example config.R value in scripts is populated with the default values used by epictope.

### License 

Epictope is distributed open-source under the GPL3 license.
