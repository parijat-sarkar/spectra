#!/usr/bin/env python3
"""
One-command launcher for Spectra.
Installs dependencies + opens browser automatically.

Usage:
    python launch.py

Then just use http://localhost:5000 in your browser.
Press Ctrl+C to stop.
"""
import subprocess
import webbrowser
import time
import sys
import os

def run():
    print("=" * 60)
    print("Spectra Launcher")
    print("=" * 60)
    print()

    print("📦 Installing dependencies (one time only)...")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "-q", "-r", "requirements.txt"],
        capture_output=True
    )

    if result.returncode != 0:
        print("❌ Failed to install dependencies:")
        print(result.stderr.decode())
        sys.exit(1)

    print("✓ Dependencies ready")
    print()
    print("🚀 Starting Spectra server...")
    print("   Opening browser in 2 seconds...")
    print()

    # Start Flask app
    server = subprocess.Popen([sys.executable, "app.py"])

    # Wait for server to start
    time.sleep(2)

    # Open browser
    try:
        webbrowser.open("http://localhost:5000")
        print("✓ Browser opened: http://localhost:5000")
    except:
        print("⚠ Could not auto-open browser. Visit: http://localhost:5000")

    print()
    print("=" * 60)
    print("Press Ctrl+C to stop the server")
    print("=" * 60)
    print()

    try:
        server.wait()
    except KeyboardInterrupt:
        print()
        print("Shutting down...")
        server.terminate()
        server.wait()
        print("✓ Done. Goodbye!")

if __name__ == "__main__":
    run()
