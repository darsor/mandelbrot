#cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)
project (MandelExplorer CXX)

# Version numbers
set (MandelExplorer_VERSION_MAJOR 0)
set (MandelExplorer_VERSION_MINOR 1)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -DUSE_SIMD_ALGORITHM")

# Adding extra libraries for displays and things
set (EXTRA_LIBS ${EXTRA_LIBS} 
    GL
    GLU
    sfml-system
    sfml-window
    sfml-graphics
)

add_executable (MandelExplorer
        mandelbrotViewer.cpp
        mandelbrotExplorer.cpp
)
target_link_libraries (MandelExplorer ${EXTRA_LIBS})
