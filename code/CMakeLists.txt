cmake_minimum_required(VERSION 3.3) # Required by SUNDIALS.
project(lotkavolterra1910)

find_package(SUNDIALS 6 REQUIRED COMPONENTS cvodes)

add_executable(lotkavolterra1910 lotkavolterra1910.c)
target_link_libraries(lotkavolterra1910 PRIVATE SUNDIALS::cvodes_shared)
