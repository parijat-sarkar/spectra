@echo off
REM Spectra launcher for Windows
REM Double-click this file to start Spectra

echo ============================================================
echo Spectra Launcher
echo ============================================================
echo.

echo Installing dependencies (one time only)...
pip install -q -r requirements.txt >nul 2>&1

if %errorlevel% neq 0 (
    echo Failed to install dependencies
    pause
    exit /b 1
)

echo Dependencies ready
echo.
echo Starting Spectra server...
echo Opening browser in 2 seconds...
echo.

REM Start the Flask app
start python app.py

REM Wait 2 seconds for server to start
timeout /t 2 /nobreak

REM Open browser
start http://localhost:5000

echo.
echo ============================================================
echo Spectra is running at http://localhost:5000
echo Close this window to stop the server
echo ============================================================
echo.

pause
