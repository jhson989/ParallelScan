# CMake bianry version
cmake_minimum_required(VERSION 3.13)



# Proejct info
project(
    ParallelBitomicSort
    VERSION 1.0
    DESCRIPTION "Parallel Bitomic Sort"
)

add_subdirectory(PrefixSum)
add_subdirectory(Scan)