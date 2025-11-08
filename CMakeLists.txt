cmake_minimum_required(VERSION 3.16)
project(mpi_battery_heat_student LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
find_package(OpenMP)

add_executable(battery_heat
    src/main.cpp
)

target_include_directories(battery_heat PRIVATE include)
target_link_libraries(battery_heat PRIVATE MPI::MPI_CXX)
if (OpenMP_CXX_FOUND)
  target_link_libraries(battery_heat PRIVATE OpenMP::OpenMP_CXX)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
