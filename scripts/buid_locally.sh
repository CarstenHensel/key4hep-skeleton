#!/bin/bash
set -e

# create install directory
mkdir -p install

# Create build directory if it doesn't exist
mkdir -p build
cd build

# Configure CMake (only needed once per clean build)
cmake .. -DCMAKE_INSTALL_PREFIX=../install -DPython_EXECUTABLE=$(which python3)

# Build everything
make -j$(nproc)

# Install to local directory
make install

echo "Local build and install completed."
