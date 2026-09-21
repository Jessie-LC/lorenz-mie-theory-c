cmakeOptions="-G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=gcc -S src -B build/linux/debug"

cmake $cmakeOptions
cmake --build build/linux/debug -v

./build/linux/debug/bin/MIE