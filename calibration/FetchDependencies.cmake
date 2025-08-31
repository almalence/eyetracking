# Fetch OpenXR SDK

set(KHR_OPENXR_SDK_REPO_DIR "${CMAKE_SOURCE_DIR}/3rd/OpenXR-SDK")

if(NOT EXISTS "${KHR_OPENXR_SDK_REPO_DIR}")
    message(STATUS "Cloning OpenXR-SDK repository...")
    execute_process(
        COMMAND git clone https://github.com/KhronosGroup/OpenXR-SDK.git ${KHR_OPENXR_SDK_REPO_DIR}
        RESULT_VARIABLE result
    )
    if(result)
        message(FATAL_ERROR "Failed to clone OpenXR-SDK repository.")
    endif()
else()
    message(STATUS "OpenXR-SDK repository already exists. Skipping clone.")
endif()

target_include_directories(EyeTrackerCalibrationApp PRIVATE 
    "${KHR_OPENXR_SDK_REPO_DIR}/src/common"
    "${KHR_OPENXR_SDK_REPO_DIR}/include/" 
    "${KHR_OPENXR_SDK_REPO_DIR}/include/openxr" 
    "${KHR_OPENXR_SDK_REPO_DIR}/src/")

set(OpenXR_SDK_BUILD_DIR "${KHR_OPENXR_SDK_REPO_DIR}/build")

# Ensure the build directory exists
if(NOT EXISTS "${OpenXR_SDK_BUILD_DIR}")
    file(MAKE_DIRECTORY "${OpenXR_SDK_BUILD_DIR}")
endif()

# Configure the OpenXR SDK
if (ANDROID)
execute_process(
    COMMAND ${CMAKE_COMMAND} -S ${KHR_OPENXR_SDK_REPO_DIR} -B ${OpenXR_SDK_BUILD_DIR} -G Ninja -DANDROID_ABI=arm64-v8a -DANDROID_PLATFORM=android-24 -DANDROID_NDK=${ANDROID_NDK} -DCMAKE_TOOLCHAIN_FILE=${ANDROID_NDK}/build/cmake/android.toolchain.cmake -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
    RESULT_VARIABLE configure_result
)
if(configure_result)
    message(FATAL_ERROR "Failed to configure OpenXR-SDK.")
endif()
else ()
    execute_process(
        COMMAND ${CMAKE_COMMAND} -S ${KHR_OPENXR_SDK_REPO_DIR} -B ${OpenXR_SDK_BUILD_DIR}
        RESULT_VARIABLE configure_result
    )
    if(configure_result)
        message(FATAL_ERROR "Failed to configure OpenXR-SDK.")
    endif()
endif ()

# Build the OpenXR SDK
execute_process(
    COMMAND ${CMAKE_COMMAND} --build ${OpenXR_SDK_BUILD_DIR}  --config ${CMAKE_BUILD_TYPE}
    RESULT_VARIABLE build_result
)
if(build_result)
    message(FATAL_ERROR "Failed to build OpenXR-SDK.")
endif()

if (ANDROID)
    target_link_libraries(EyeTrackerCalibrationApp PRIVATE 
        "${OpenXR_SDK_BUILD_DIR}/src/loader/libopenxr_loader.so"
    )
else()
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(OPENXR_LOADER_LIB_NAME "openxr_loaderd.lib")
    else()
        set(OPENXR_LOADER_LIB_NAME "openxr_loader.lib")
    endif()
    target_link_libraries(EyeTrackerCalibrationApp PRIVATE 
        "${OpenXR_SDK_BUILD_DIR}/src/loader/${CMAKE_BUILD_TYPE}/${OPENXR_LOADER_LIB_NAME}"
    )
endif()

# Fetch some OpenXR sources

set(OPENGL_DEPENDENCIES_DIR "${CMAKE_SOURCE_DIR}/3rd/OpenGL")
if(NOT EXISTS "${OPENGL_DEPENDENCIES_DIR}")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/refs/heads/main/src/common/gfxwrapper_opengl.h" 
        "${OPENGL_DEPENDENCIES_DIR}/gfxwrapper_opengl.h")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/refs/heads/main/src/common/gfxwrapper_opengl.c" 
        "${OPENGL_DEPENDENCIES_DIR}/gfxwrapper_opengl.c")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/external/include/GL/glext.h" 
        "${OPENGL_DEPENDENCIES_DIR}/GL/glext.h")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/external/include/GL/glext.h.license"
            "${OPENGL_DEPENDENCIES_DIR}/GL/glext.h.license")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/external/include/GL/gl_format.h" 
        "${OPENGL_DEPENDENCIES_DIR}/GL/gl_format.h")

    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/src/external/glad2/include/EGL/eglplatform.h"
        "${OPENGL_DEPENDENCIES_DIR}/EGL/eglplatform.h")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/src/external/glad2/include/KHR/khrplatform.h"
        "${OPENGL_DEPENDENCIES_DIR}/KHR/khrplatform.h")

    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/src/external/glad2/include/glad/gl.h"
        "${OPENGL_DEPENDENCIES_DIR}/glad/gl.h")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/src/external/glad2/src/gl.c"
        "${OPENGL_DEPENDENCIES_DIR}/gl.c")

    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/src/external/glad2/include/glad/egl.h"
        "${OPENGL_DEPENDENCIES_DIR}/glad/egl.h")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/src/external/glad2/src/egl.c"
        "${OPENGL_DEPENDENCIES_DIR}/egl.c")

    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/src/external/glad2/include/glad/wgl.h"
        "${OPENGL_DEPENDENCIES_DIR}/glad/wgl.h")
    file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/src/external/glad2/src/wgl.c"
        "${OPENGL_DEPENDENCIES_DIR}/wgl.c")

    if (WIN32)
        file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/external/include/GL/wglext.h"
                "${OPENGL_DEPENDENCIES_DIR}/GL/wglext.h")
        file(DOWNLOAD "https://raw.githubusercontent.com/KhronosGroup/OpenXR-SDK-Source/main/external/include/GL/wglext.h.license"
            "${OPENGL_DEPENDENCIES_DIR}/GL/wglext.h.license")
    endif ()
else()
    message(STATUS "Requred OpenGL files already exists. Skipping downloading.")
endif()

add_library(openxr-gfxwrapper STATIC ${OPENGL_DEPENDENCIES_DIR}/gfxwrapper_opengl.c ${OPENGL_DEPENDENCIES_DIR}/gl.c ${OPENGL_DEPENDENCIES_DIR}/wgl.c)
target_include_directories(openxr-gfxwrapper PRIVATE "${OPENGL_DEPENDENCIES_DIR}")
target_include_directories(EyeTrackerCalibrationApp PRIVATE "${OPENGL_DEPENDENCIES_DIR}")
target_link_libraries(EyeTrackerCalibrationApp PRIVATE openxr-gfxwrapper)
