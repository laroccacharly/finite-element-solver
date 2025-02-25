#!/bin/bash

# Exit on error
set -e

# Detect Mac architecture and set Qt path accordingly
if [[ $(uname -m) == 'arm64' ]]; then
    # Apple Silicon (M1/M2/M3)
    QT_PATH="/opt/homebrew/opt/qt@5"
else
    # Intel Mac
    QT_PATH="/usr/local/opt/qt@5"
fi

echo "Using Qt path: $QT_PATH"

# Create build directory
mkdir -p build
cd build

# Configure with CMake
cmake -DCMAKE_PREFIX_PATH=$QT_PATH -DCMAKE_BUILD_TYPE=Release ..

# Build
make -j4

# Run the application
./finite-element-solver 