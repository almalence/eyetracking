# Desktop-specific configurations

# Add executable for desktop build
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

find_package(OpenGL REQUIRED)

find_package(OpenCV REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS} )

target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL ${OpenCV_LIBS} )

