cmake_minimum_required(VERSION 3.8)
project(AudioEdit VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)

add_executable(audio_edit src/main.cpp)