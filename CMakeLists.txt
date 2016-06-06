project(vad)
cmake_minimum_required(VERSION 2.8)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMAKE_HELPERS/")

set(HEADER_FILES ${PROJECT_SOURCE_DIR}/vad.h)

aux_source_directory(. SRC_LIST)
include_directories(${YOUR_DIRECTORY})

find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)
find_package(FFTW REQUIRED)
include_directories(FFTW_INCLUDE_DIR)

add_executable(vad voice_detection.cpp )
target_link_libraries(vad -lm -lfftw3 -lstdc++ )
