@echo off

rd /s /q "build/windows/debug"
set cmakeOptions=-G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=gcc -S src -B build/windows/debug

cmake %cmakeOptions%
cmake --build build/windows/debug

"./build/windows/debug/bin/MIE.exe"
