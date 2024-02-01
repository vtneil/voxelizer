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

add_executable(voxelizer_bvh
        ${SOURCES}
        main_bvh.cpp
)

add_executable(voxelizer_bvh_fp
        ${SOURCES}
        main_bvh.cpp
)

target_compile_options(voxelizer PRIVATE
        $<$<CONFIG:Release>:-O3 -mavx512f -mavx512cd -mavx512er -mavx512pf>
)

target_compile_options(voxelizer_bvh PRIVATE
        $<$<CONFIG:Release>:-O3 -mavx512f -mavx512cd -mavx512er -mavx512pf>
)

target_compile_options(voxelizer_bvh_fp PRIVATE
        $<$<CONFIG:Release>:-O3 -mavx512f -mavx512cd -mavx512er -mavx512pf>
)
