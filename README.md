# Epitope Tag

Simplified code for identifying ideal epitope tag insertion sites for proteins.

## Table of Contents

- [Config](#config)
- [Installation](#installation)
- [Dependencies](#dependencies)

## Config

`install.R` 

## Installation

A installation script installs necessary R packages, downloads CDS files from NCBI FTP, and creates protein blastable databases.

```
# Rscript install.R
```
## Dependencies

R packages
- UniProt.ws
- httr
- rMSA
- Biostrings
- rvest

Programs
- blast-plus/2.13.0
- muscle/3.8.1551
- dssp/3.1.4
- r/4.2.2

## Usage

Set the epitope tag insertion site parameters in provided `config` file. 

```
# Rscript single_score.R
```

Folder Structure
- `code/`: 
- `data/`: 
  - `models/`: Folder for downloaded query pdb files
  - `CDS`: Folder for CDS files for species in MSA
- `outputs/`: 
