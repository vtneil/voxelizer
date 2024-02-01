cmake_minimum_required(VERSION 3.22)
project(voxelizer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 17)

include_directories(include)
file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(voxelizer
        ${SOURCES}
        main.cpp
)

add_executable(voxelizer_vdb
        ${SOURCES}
        main_vdb.cpp
)

add_executable(voxelizer_test
        ${SOURCES}
        test.cpp
)

target_compile_options(voxelizer PRIVATE
        -O3
)
