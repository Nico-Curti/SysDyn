#!/bin/bash

if [[ "$OSTYPE" == "darwin"* ]]; then
 export CC="/usr/local/bin/gcc-7"
 export CXX="/usr/local/bin/g++-7"
fi

rm -rf build
mkdir -p build
cd build

#sudo apt-get install ninja-build
cmake -G "Ninja" "-DCMAKE_BUILD_TYPE=Release" ..
cmake --build .
cd ..

