#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include "vsf.h"

using namespace std;

PYBIND11_MODULE(stats_library, m) {
    m.doc() = string("Module that regroups the necessary statistic and analysis tools to compute n-th order structure"
                     "functions");
    m.def("str_func_cpp", &structure_function,
          "Compute the n-th order structure function of a two-dimensional array.");
}

// MAC :    clang++ -std=c++17 -shared -undefined dynamic_lookup -I./pybind11/include/ `python3.12 -m pybind11 --includes` vsf.cpp stats.cpp tools.cpp time.cpp pybind11.cpp -o stats_library.so `python3.12-config --ldflags` -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp
// LINUX :  g++ -std=c++17 -shared -fPIC -I./pybind11/include/ `python3.12 -m pybind11 --includes` vsf.cpp stats.cpp tools.cpp time.cpp pybind11.cpp -o stats_library.so `python3.12-config --ldflags` -fopenmp -lm
// C++:     clang++ -std=c++17 vsf.cpp stats.cpp tools.cpp time.cpp -o test -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp
