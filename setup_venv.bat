@echo off
rem Create the Python environment for the toolkit and install the requirements.
rem
rem   setup_venv.bat        Anaconda / Miniconda found: conda environment in .conda_env
rem                         (Python 3.12 from conda-forge with a current C++ runtime,
rem                         so PySide6 loads), otherwise a plain venv in .venv
rem   setup_venv.bat venv   always use a plain venv in .venv
rem
rem Run it again to update or repair the environment.
setlocal
cd /d "%~dp0"
set "ENVDIR=%~dp0.conda_env"
set "VPY="
if /i "%~1"=="venv" goto :venv

rem ---------------------------------------------------------------- conda --
if exist "%ENVDIR%\python.exe" goto :conda_ready
call :find_conda
if not defined CONDA goto :venv
for %%C in ("%CONDA%") do set "CROOT=%%~dpC"
for %%C in ("%CROOT%..") do set "CROOT=%%~fC"
if exist "%CROOT%\condabin\conda.bat" set "CONDA=%CROOT%\condabin\conda.bat"
echo Using conda: "%CONDA%"
echo Creating conda environment .conda_env (Python 3.12 from conda-forge) ...
setlocal
set "PATH=%CROOT%;%CROOT%\Library\bin;%CROOT%\Scripts;%CROOT%\condabin;%PATH%"
call "%CONDA%" create -y -p "%ENVDIR%" --override-channels -c conda-forge python=3.12 pip
endlocal
if not exist "%ENVDIR%\python.exe" (
    echo conda could not create the environment - using a plain venv instead.
    goto :venv
)
:conda_ready
rem minimal activation: the environment's own DLLs come first
set "PATH=%ENVDIR%;%ENVDIR%\Library\bin;%ENVDIR%\Scripts;%PATH%"
set "VPY=%ENVDIR%\python.exe"
goto :install

rem ----------------------------------------------------------------- venv --
:venv
set "PY=python"
where py >nul 2>nul && (
    set "PY=py -3"
    for %%V in (3.13 3.11 3.12) do py -%%V -c "" >nul 2>nul && set "PY=py -%%V"
)
if not exist .venv\Scripts\python.exe (
    echo Creating virtual environment .venv with: %PY%
    %PY% -m venv .venv || goto :fail
)
set "VPY=.venv\Scripts\python.exe"

rem -------------------------------------------------------------- install --
:install
"%VPY%" -c "import sys; print('Python', sys.version.split()[0], sys.executable)"
"%VPY%" -m pip install --upgrade pip
"%VPY%" -m pip install -r requirements.txt || goto :fail

rem check that Qt (PySide6) can be loaded
"%VPY%" run_gui.py --check && goto :ok
if /i not "%VPY%"==".venv\Scripts\python.exe" goto :fail
rem venv only: older PySide6 builds (< 6.8) accept older C++ runtimes, Python <= 3.12
"%VPY%" -c "import sys; sys.exit(0 if sys.version_info < (3, 13) else 1)" || goto :newpy
echo.
echo Trying an older PySide6 that works with older C++ runtimes ...
"%VPY%" -m pip install "PySide6<6.8" || goto :fail
"%VPY%" run_gui.py --check && goto :ok
goto :fail

:newpy
echo.
echo Qt cannot be loaded, see the diagnosis above. With this Python version
echo no older PySide6 exists. Install Anaconda/Miniconda or Python 3.12 and
echo run setup_venv.bat again.
goto :fail

:ok
echo.
echo Done. Start the GUI with run_gui.bat
exit /b 0
:fail
echo Setup failed - see the messages above.
pause
exit /b 1

rem ------------------------------------------------------------------------
:find_conda
rem conda.exe / conda.bat: active conda shell, PATH, python on PATH, usual folders
set "CONDA="
if defined CONDA_EXE (
    if exist "%CONDA_EXE%" set "CONDA=%CONDA_EXE%"
)
if defined CONDA exit /b 0
for /f "delims=" %%P in ('where conda 2^>nul') do (
    if not defined CONDA set "CONDA=%%P"
)
if defined CONDA exit /b 0
for /f "delims=" %%P in ('where python 2^>nul') do (
    if not defined CONDA if exist "%%~dpPScripts\conda.exe" set "CONDA=%%~dpPScripts\conda.exe"
)
if defined CONDA exit /b 0
for %%R in ("%USERPROFILE%\anaconda3" "%USERPROFILE%\miniconda3" "%LOCALAPPDATA%\anaconda3" "%LOCALAPPDATA%\miniconda3" "%ProgramData%\anaconda3" "%ProgramData%\miniconda3") do (
    if not defined CONDA if exist "%%~R\Scripts\conda.exe" set "CONDA=%%~R\Scripts\conda.exe"
)
exit /b 0
