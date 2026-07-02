cmake_minimum_required(VERSION 3.25.0)
project(EyeTracking VERSION 0.1.0 LANGUAGES C)

add_compile_options(-DHAVE_AVX2 -ffast-math -ftree-vectorize -mavx2 -mfpmath=sse)

# static library
add_library(EyeTrack EyeTracking.c EyeCalibration.c ellipse.c eigen.c)

# cmd-line tracking test
add_executable(EyeTracking EyeTracking.c ellipse.c eigen.c test.c)
if (WIN32)
else()
target_link_libraries (EyeTracking m)
endif()

# cmd-line calibration
add_executable(EyeCalibration EyeTracking.c EyeCalibration.c ellipse.c eigen.c calibration.c)
if (WIN32)
else()
target_link_libraries (EyeCalibration m)
endif()

