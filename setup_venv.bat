@echo off
rem Create .venv next to this file and install the requirements.
rem Usage: setup_venv.bat   (run again to update)
setlocal
cd /d "%~dp0"
set "PY=python"
where py >nul 2>nul && set "PY=py -3"
if not exist .venv\Scripts\python.exe (
    echo Creating virtual environment .venv ...
    %PY% -m venv .venv || goto :fail
)
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -r requirements.txt || goto :fail
echo.
echo Done. Start the GUI with run_gui.bat
exit /b 0
:fail
echo Setup failed - see the messages above.
exit /b 1
