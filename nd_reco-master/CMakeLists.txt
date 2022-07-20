cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B1)
set(CMAKE_BUILD_TYPE Debug)
find_package(EDepSim)
find_package(ROOT REQUIRED COMPONENTS EG)
include("${ROOT_USE_FILE}")

include_directories(${PROJECT_SOURCE_DIR})

add_executable(smearBMgo smearBeamMonitoring_SAND.C)
target_link_libraries(smearBMgo PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})


