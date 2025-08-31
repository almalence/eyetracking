# Android-specific configurations

# Add library for Android build (usually shared library for Android NDK)
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})

# Android specific configurations
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DANDROID")

# Link libraries for Android
include_directories(${ANDROID_NDK}/sources/android/native_app_glue)
add_library(android_native_app_glue STATIC ${ANDROID_NDK}/sources/android/native_app_glue/android_native_app_glue.c)

target_link_libraries(${PROJECT_NAME} PUBLIC android_native_app_glue log android)

# Find and link OpenGL ES libraries based on your requirements
find_library(GLESv3_LIB GLESv3) # Use GLESv3 if targeting OpenGL ES 3.0
find_library(EGL_LIB EGL)

find_package(OpenCV REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS} )

# Link OpenGL ES, EGL, OpenCV libraries to the main project target
target_link_libraries(${PROJECT_NAME} PRIVATE ${GLESv3_LIB} ${EGL_LIB} ${OpenCV_LIBS} )
target_link_options(${PROJECT_NAME} PRIVATE "-uANativeActivity_onCreate")
