@echo off
conda config --set always_yes yes
conda activate epictope
conda install -c bioconda blast muscle
conda install -c speleo3 dssp
conda install -c conda-forge r-base r-stringi r-openssl
conda install -c conda-forge r-base r-stringi r-openssl gzip git
conda config --set always_yes no
