#!/bin/bash
conda config --set always_yes yes
conda activate epictope
conda install -c bioconda blast
conda install -c bioconda muscle
conda install -c salilab dssp
conda install -c conda-forge r-base
conda install -c conda-forge r-stringi
conda install -c conda-forge r-openssl
conda config --set always_yes no
