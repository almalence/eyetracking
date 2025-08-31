cmake_minimum_required(VERSION 3.11)

project(EyeTrackerCalibrationApp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCE_FILES
    src/common.h
    src/main.cpp
    src/openxr_module.h
    src/openxr_module.cpp
    src/graphics_module.h
    src/graphics_module.cpp
    src/shaders.h
    src/interaction_module.h
    src/interaction_module.cpp
    src/capture_module.h
    src/capture_module.cpp
)

################################# TARGETS #################################

# Include platform-specific configuration
if (ANDROID)
    add_compile_definitions(XR_USE_PLATFORM_ANDROID)
    add_compile_definitions(XR_USE_PLATFORM_EGL)
    add_compile_definitions(XR_USE_GRAPHICS_API_OPENGL_ES)
    include(${CMAKE_SOURCE_DIR}/Android.cmake)
elseif (WIN32)
    add_compile_definitions(XR_USE_PLATFORM_WIN32)
    add_compile_definitions(XR_USE_GRAPHICS_API_OPENGL)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
    include(${CMAKE_SOURCE_DIR}/Desktop.cmake)
endif()

include(${CMAKE_SOURCE_DIR}/FetchDependencies.cmake)

# Eye Tracking lib
target_include_directories(${PROJECT_NAME} PRIVATE et-lib)
target_link_directories(${PROJECT_NAME} PRIVATE et-lib)
target_link_libraries(${PROJECT_NAME} PRIVATE "libEyeTrack")


############################# CMAKE REFERENCES #############################
set(CMAKE_FILES
    Android.cmake
    Desktop.cmake
    FetchDependencies.cmake
)

# Create a custom target to include the files in the solution
add_custom_target(CMakeProjectFiles ALL SOURCES ${CMAKE_FILES})

# Include the custom target in the main solution
add_dependencies(EyeTrackerCalibrationApp CMakeProjectFiles)
