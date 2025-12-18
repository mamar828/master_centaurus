#!/bin/bash

# Build script for stats_library pybind11 module

# Create build directory if it doesn't exist
mkdir -p build
cd build

# Configure with CMake
cmake ..

# Build
cmake --build .

echo "Build complete! stats_library.so is ready."
