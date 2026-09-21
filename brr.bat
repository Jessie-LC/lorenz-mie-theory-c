@echo off

rd /s /q "build/windows/release"
set cmakeOptions=-G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang -S src -B build/windows/release

cmake %cmakeOptions%
cmake --build build/windows/release

"./build/windows/release/bin/MIE.exe"
