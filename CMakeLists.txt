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

add_executable(voxelizer_cpu
        ${SOURCES}
        main_bvh_omp_avx.cpp
)

target_compile_options(voxelizer PRIVATE
        -O3
)

target_compile_options(voxelizer_bvh PRIVATE
        -O3
)

target_compile_options(voxelizer_cpu PRIVATE
        -O3 -fopenmp -flto -mavx512f -mavx512cd -mavx512er -mavx512pf
)
target_link_libraries(voxelizer_cpu PRIVATE
        -O3 -fopenmp -flto -mavx512f -mavx512cd -mavx512er -mavx512pf
)

target_compile_options(voxelizer_bvh_fp PRIVATE
        -O3
)
