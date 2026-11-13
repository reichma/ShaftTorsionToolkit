@echo off
rem Create .venv next to this file and install the requirements.
rem Usage: setup_venv.bat   (run again to update or repair)
setlocal
cd /d "%~dp0"
set "PY=python"
where py >nul 2>nul && set "PY=py -3"
if not exist .venv\Scripts\python.exe (
    echo Creating virtual environment .venv ...
    %PY% -m venv .venv || goto :fail
)
set "VPY=.venv\Scripts\python.exe"
%VPY% -m pip install --upgrade pip
%VPY% -m pip install -r requirements.txt || goto :fail

rem --- check that Qt (PySide6) can be loaded --------------------------------
%VPY% -c "import PySide6.QtCore" >nul 2>nul && goto :ok
echo.
echo PySide6 is installed but Qt cannot be loaded (DLL load failed).
echo This usually means the Python behind .venv (e.g. Anaconda) ships an
echo older C++ runtime than the newest PySide6 needs.
echo Installing an older PySide6 that works with it ...
%VPY% -m pip install "PySide6<6.8" || goto :fail
%VPY% -c "import PySide6.QtCore" && goto :ok
echo.
echo Qt still cannot be loaded. Run run_gui.bat to see a detailed diagnosis.
goto :fail

:ok
echo.
echo Done. Start the GUI with run_gui.bat
exit /b 0
:fail
echo Setup failed - see the messages above.
pause
exit /b 1
