cmake_minimum_required(VERSION 3.12)
project(a_birch LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(birch_corrected
    main.cpp
    Common.cpp
    Dataset.cpp
    KNNOutlier.cpp
    BIRCH.cpp
    Validation.cpp
    Pipeline.cpp)

if(MSVC)
    target_compile_options(birch_corrected PRIVATE /W4)
else()
    target_compile_options(birch_corrected PRIVATE -O2 -Wall -Wextra -Wpedantic)
endif()
