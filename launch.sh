#!/bin/bash
# Spectra launcher for Mac/Linux
# Double-click or run: bash launch.sh

echo "============================================================"
echo "Spectra Launcher"
echo "============================================================"
echo ""

echo "📦 Installing dependencies (one time only)..."
pip install -q -r requirements.txt

if [ $? -ne 0 ]; then
    echo "❌ Failed to install dependencies"
    exit 1
fi

echo "✓ Dependencies ready"
echo ""
echo "🚀 Starting Spectra server..."
echo "   Opening browser in 2 seconds..."
echo ""

# Start the Flask app
python app.py &
FLASK_PID=$!

# Wait for server to start
sleep 2

# Open browser
open "http://localhost:5000" 2>/dev/null || xdg-open "http://localhost:5000" 2>/dev/null

echo "✓ Browser opened: http://localhost:5000"
echo ""
echo "============================================================"
echo "Press Ctrl+C to stop the server"
echo "============================================================"
echo ""

# Wait for Flask process
wait $FLASK_PID
