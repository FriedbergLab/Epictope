#!/bin/bash
# Install conda
conda config --set always_yes yes
conda create -n epictope
conda activate epictopeR
conda install -c bioconda blast muscle
conda install -c salilab dssp
conda install -c conda-forge r-base r-stringi r-openssl r-remotes
conda config --set always_yes no

# Install R packages
R -e "remotes::install_github('henrichung/epitope_tag')"
# Install epitope_tag scripts
wget "https://raw.githubusercontent.com/henrichung/epitope_tag/main/scripts/single_score.R" 
wget "https://raw.githubusercontent.com/henrichung/epitope_tag/main/scripts/plot_scores.R" 
wget "https://raw.githubusercontent.com/henrichung/epitope_tag/main/scripts/install.R"
wget "https://raw.githubusercontent.com/henrichung/epitope_tag/main/scripts/config_defaults.R" 