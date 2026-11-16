@echo off
rem Start the Shaft Torsion Toolkit GUI. Optional argument: a .step file or a project .json
cd /d "%~dp0"
if exist .conda_env\python.exe goto :run
if exist .venv\Scripts\python.exe goto :run
call setup_venv.bat || exit /b 1
:run
if exist .conda_env\python.exe (
    set "PATH=%~dp0.conda_env;%~dp0.conda_env\Library\bin;%~dp0.conda_env\Scripts;%PATH%"
    set "VPY=.conda_env\python.exe"
) else (
    set "VPY=.venv\Scripts\python.exe"
)
"%VPY%" run_gui.py %*
if errorlevel 1 (
    echo.
    echo The GUI stopped with an error - see the messages above.
    pause
)
