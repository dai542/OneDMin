#!/usr/bin/env bash

mkdir -p build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_Fortran_COMPILER=$FC

make VERBOSE=1
make install
