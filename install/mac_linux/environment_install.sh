#!/bin/bash
conda config --set always_yes yes
conda create -n epictope
conda activate epictope
conda install -c bioconda blast muscle
conda install -c salilab dssp
conda install -c conda-forge r-base r-stringi r-openssl r-remotes git
conda config --set always_yes no
