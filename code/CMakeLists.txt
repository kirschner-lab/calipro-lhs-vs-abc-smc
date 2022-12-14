cmake_minimum_required(VERSION 3.3) # Required by SUNDIALS.
project(loktavolterra1910)

find_package(SUNDIALS 6 REQUIRED COMPONENTS cvodes)

add_executable(loktavolterra1910 loktavolterra1910.c)
target_link_libraries(loktavolterra1910 PRIVATE SUNDIALS::cvodes_shared)
