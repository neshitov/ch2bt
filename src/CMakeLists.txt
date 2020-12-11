cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(snf C)

# add the executable


find_package(OpenMP REQUIRED)
add_executable(gap_snf gap_snf.c snf.c)
target_link_libraries(gap_snf PUBLIC  OpenMP::OpenMP)
target_compile_options(gap_snf PUBLIC OpenMP::OpenMP)
