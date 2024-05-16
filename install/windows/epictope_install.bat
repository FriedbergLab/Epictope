@echo off
:: CONDA INSTALLATION
:: check if conda is installed
where conda >nul 2>nul
if %errorlevel% neq 0 (
    echo Conda is not installed on your system.
    echo Please install it first.
    exit /b
)
:: create the 'epictope' environment if it doesn't exist
call conda info --envs | findstr /r "^epictope " >nul
if %errorlevel% neq 0 (
    call conda create --name epictope --yes
)
:: activate the 'Epictope' environment and install dependencies
call conda activate epictope
call conda install -c speleo3 dssp
call conda install -c conda-forge r-base r-stringi r-openssl r-remotes
:: report success
echo Epictope environment installation complete.

:: BLAST INSTALLATION
:: Check if the file exists
if exist ncbi-blast-2.7.1+-win64.exe (
    echo Blast executable found. Skipping installation.
) else (
    :: Download the blast installer
    echo Downloading BLAST+ installer...
    powershell -Command "(New-Object Net.WebClient).DownloadFile('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-win64.exe', 'ncbi-blast-2.7.1+-win64.exe')"
    :: Run the blast installer
    echo Running BLAST+ installer...
    start /wait ncbi-blast-2.7.1+-win64.exe /S
    :: Add blast to path
    echo Adding BLAST+ to the PATH...
    set "PATH=%PATH%;C:\Program Files\NCBI\blast-2.7.1+\bin"
    :: report success
    echo BLAST+ installation complete.
)

:: MUSCLE INSTALLATION
:: Check if the file exists
if exist muscle.exe (
    echo MUSCLE executable found. Skipping installation.
) else (
    :: Download the MUSCLE installer
    echo Downloading MUSCLE...
    powershell -Command "(New-Object Net.WebClient).DownloadFile('https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.win64.exe', '%cd%\muscle5.1.win64.exe')"
    echo Renaming MUSCLE executable...
    move muscle5.1.win64.exe muscle.exe
    :: report success
    echo MUSCLE installation complete.
)

:: EPICTOPE INSTALLATION
:: Check if the file exists
call conda activate epictope
call R -e "remotes::install_github('FriedbergLab/EpicTope')"
call curl -o "single_score.R" "https://raw.githubusercontent.com/FriedbergLab/EpicTope/main/scripts/single_score.R"
call curl -o "plot_scores.R" "https://github.com/FriedbergLab/EpicTope/blob/main/scripts/plot_scores.R"
call curl -o "install.R" "https://raw.githubusercontent.com/FriedbergLab/EpicTope/main/scripts/install.R"
call curl -o "config_defaults.R" "https://raw.githubusercontent.com/FriedbergLab/EpicTope/main/scripts/config_defaults.R"
