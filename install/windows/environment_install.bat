@echo off

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

:: activate the 'Epictope' environment
call conda activate epictope

:: install dependencies
call conda install -c speleo3 dssp
call conda install -c conda-forge r-base r-stringi r-openssl r-remotes git

:: report success
echo.
echo Epictope environment installation complete.
exit /b
