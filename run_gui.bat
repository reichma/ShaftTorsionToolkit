@echo off
rem Start the Shaft Torsion Toolkit GUI. Optional argument: a .step file or a project .json
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe call setup_venv.bat
.venv\Scripts\python.exe run_gui.py %*
