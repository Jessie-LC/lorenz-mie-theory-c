cmake_minimum_required(VERSION 3.30)
project(LORENZMIETHEORY)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_compile_definitions(UNICODE _UNICODE _CRT_SECURE_NO_WARNINGS)

include_directories(
    AFTER
    "${CMAKE_SOURCE_DIR}/../external/khronos"
)

add_executable(MIE)
target_compile_definitions(MIE PUBLIC c_std_23)


if(WIN32)


if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(MIE PRIVATE -mavx2 -mfma -O3 -fno-fast-math)
    target_link_libraries(
        MIE PUBLIC
        opengl32.lib
    )
else(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(MIE PRIVATE -O0 -mavx2)
    target_link_libraries(
        MIE PUBLIC
        opengl32.lib
    )
endif()

target_link_directories(
    MIE
    PRIVATE
    "C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0/ucrt/x64"
    "C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0/um/x64"
)

target_include_directories(
    MIE
    PRIVATE
    "/usr/share/libc++/v1"
    "C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0/ucrt/x64"
    "C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0/um/x64"
    "../include"
)


else()


if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(MIE PRIVATE -mavx2 -mfma -O3 -fno-fast-math)
    target_link_libraries(
        MIE PUBLIC
    )
else(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(MIE PRIVATE -O0 -mavx2)
    target_link_libraries(
        MIE PUBLIC
    )
endif()

target_link_directories(
    MIE
    PRIVATE
)

target_include_directories(
    MIE
    PRIVATE
    "/usr/share/libc++/v1"
    "../include"
)

endif()

target_sources(
    MIE
    PRIVATE
    main.c
)
