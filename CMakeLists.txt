cmake_minimum_required(VERSION 3.16)
project(FRNS_bootstrap)

set(CMAKE_CXX_STANDARD 17)
file(GLOB SRC "./src/*.cpp")
file(GLOB INCLUDE "./src/*.h")
file(GLOB DATA "./data/*.h")
file(GLOB RUN "./run/*.cpp")

include_directories("./src/")

foreach (app ${RUN})
    get_filename_component(exe ${app} NAME_WE)
    add_executable(${exe} ${app} ${SRC} ${INCLUDE} ${DATA})
endforeach ()