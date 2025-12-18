cmake_minimum_required(VERSION 3.15)
project(stats_library)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find pybind11 - try multiple methods
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    # Try to use pybind11 from Python package
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
        OUTPUT_VARIABLE pybind11_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE pybind11_RESULT
    )
    if(pybind11_RESULT EQUAL 0)
        find_package(pybind11 CONFIG REQUIRED)
    else()
        message(FATAL_ERROR "pybind11 not found. Install with: pip install pybind11")
    endif()
endif()

# Find OpenMP with hints for macOS Homebrew
if(APPLE)
    # Try to find OpenMP in Homebrew location
    set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()

find_package(OpenMP REQUIRED)

# Create the Python module
pybind11_add_module(stats_library
    pybind11.cpp
    vsf.cpp
    stats.cpp
    tools.cpp
)

# Link OpenMP
target_link_libraries(stats_library PRIVATE OpenMP::OpenMP_CXX)

# Set output name (removes the default prefix/suffix that might be added)
set_target_properties(stats_library PROPERTIES
    PREFIX ""
    OUTPUT_NAME "stats_library"
)

# Platform-specific settings
if(APPLE)
    set_target_properties(stats_library PROPERTIES
        SUFFIX ".so"
        LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()
