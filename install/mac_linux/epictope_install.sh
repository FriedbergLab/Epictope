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

    # Create a new conda environment for epictope
    conda create -n epictope
    conda install -n epictope -c bioconda blast muscle
    conda install -n epictope -c salilab dssp
    conda install -n epictope -c anaconda "openssl>=3.0.9"
    conda install -n epictope -c conda-forge r-base r-stringi r-openssl r-remotes "python>=3.11.4"
    conda install -n epictope libboost=1.73.0 # for compatibility for dssp 3
    
    # Install R packages in the epictope environment
    (
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate epictope
    R -e "remotes::install_github('FriedbergLab/EpicTope')"
    conda deactivate
    )
    # Install epitope_tag scripts
    curl -O "https://raw.githubusercontent.com/FriedbergLab/Epictope/main/scripts/single_score.R" 
    curl -O "https://raw.githubusercontent.com/FriedbergLab/Epictope/main/scripts/plot_scores.R" 
    curl -O "https://raw.githubusercontent.com/FriedbergLab/Epictope/main/scripts/install.R"
    curl -O "https://raw.githubusercontent.com/FriedbergLab/Epictope/main/scripts/config_defaults.R" 
fi



