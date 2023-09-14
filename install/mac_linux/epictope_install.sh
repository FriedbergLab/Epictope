#!/bin/bash

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed, installing..."

    # Determine the OS and set the appropriate miniconda installer URL
    case "$(uname -s)" in
        Darwin)
            installer_url="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
            installer_name="Miniconda3-latest-MacOSX-x86_64.sh"
            ;;
        Linux)
            installer_url="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
            installer_name="Miniconda3-latest-Linux-x86_64.sh"
            ;;
        *)
            echo "Unsupported operating system."
            exit 1
            ;;
    esac
    # Download and install conda 
    curl -O $installer_url
    bash $installer_name
    rm $installer_name
    # Initialize 
    source ~/.bashrc
    echo "Conda installed successfully. Please close and restart your terminal."
else
    echo "Conda is already installed."
fi


# Create a new conda environment for epictope
conda create -n epictope
conda activate epictope
conda install -c bioconda blast muscle
conda install -c salilab dssp
conda install -c anaconda openssl>=3.0.9
conda install -c conda-forge r-base r-stringi r-openssl r-remotes python>=3.11.4

# Install R packages
R -e "remotes::install_github('henrichung/epitope_tag')"
# Install epitope_tag scripts
wget "https://raw.githubusercontent.com/henrichung/epitope_tag/main/scripts/single_score.R" 
wget "https://raw.githubusercontent.com/henrichung/epitope_tag/main/scripts/plot_scores.R" 
wget "https://raw.githubusercontent.com/henrichung/epitope_tag/main/scripts/install.R"
wget "https://raw.githubusercontent.com/henrichung/epitope_tag/main/scripts/config_defaults.R" 
