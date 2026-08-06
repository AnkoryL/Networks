@echo off

title Networks - Gene Network Analysis
echo.
echo ==========================================
echo          NETWORKS
echo ==========================================
echo.
echo Starting application.
echo.
echo Loading R and libraries. Please wait.
echo.

set "ROOT=%~dp0"
set "APP=%ROOT%app"
set "R_HOME=%APP%\portable-r-4.6.1-win-x64"

cd /d "%APP%"

"%R_HOME%\bin\Rscript.exe" "%APP%\app.R"


