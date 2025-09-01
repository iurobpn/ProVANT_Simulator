#!/bin/bash

wget https://cmake.org/files/v3.23/cmake-3.23.5.tar.gz
tar -xvf cmake-3.23.5.tar.gz
cd cmake-3.23.5

./bootstrap --system-curl
make -j$(nproc)
make install

cd ..
rm -rf cmake-3.23.5

