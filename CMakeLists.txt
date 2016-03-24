cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)
project (MandelExplorer CXX)

# Version numbers
set (MandelExplorer_VERSION_MAJOR 0)
set (MandelExplorer_VERSION_MINOR 1)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O3")

# Adding extra libraries for displays and things
set (EXTRA_LIBS ${EXTRA_LIBS} 
    GL
    GLU
    mpfr
    sfml-system
    sfml-window
    sfml-graphics
    pthread
)

add_executable (MandelExplorer
        mandelbrotViewer.cpp
        mandelbrotExplorer.cpp
)
target_link_libraries (MandelExplorer ${EXTRA_LIBS})
