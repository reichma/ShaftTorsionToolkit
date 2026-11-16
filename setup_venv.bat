@echo off
rem Create .venv next to this file and install the requirements.
rem Usage: setup_venv.bat   (run again to update or repair)
setlocal
cd /d "%~dp0"

rem --- pick the Python: prefer 3.12 / 3.11 / 3.13 from the py launcher -------
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
%VPY% -c "import sys; print('Python', sys.version.split()[0], sys.executable)"
%VPY% -m pip install --upgrade pip
%VPY% -m pip install -r requirements.txt || goto :fail

rem --- check that Qt (PySide6) can be loaded --------------------------------
%VPY% run_gui.py --check && goto :ok
rem Older PySide6 builds (< 6.8) need a less recent C++ runtime but exist only
rem for Python <= 3.12.
%VPY% -c "import sys; sys.exit(0 if sys.version_info < (3, 13) else 1)" || goto :newpy
echo.
echo Trying an older PySide6 that works with older C++ runtimes ...
%VPY% -m pip install "PySide6<6.8" || goto :fail
%VPY% run_gui.py --check && goto :ok
goto :fail

:newpy
echo.
echo Qt cannot be loaded, see the diagnosis above. With this Python version
echo no older PySide6 exists, so the automatic fallback is not possible.
echo Recommended: install Python 3.12 (64-bit) from python.org, delete the
echo .venv folder and run setup_venv.bat again - it prefers Python 3.12.
goto :fail

:ok
echo.
echo Done. Start the GUI with run_gui.bat
exit /b 0
:fail
echo Setup failed - see the messages above.
pause
exit /b 1
