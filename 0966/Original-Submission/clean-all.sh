#!/bin/bash

cd art-gallery-pg
./clean-cmake.sh
cd -

cd grid
./clean-cmake.sh
cd -

cd polygon
./clean-cmake.sh
cd -

cd scp-solver
./clean-cmake.sh
cd -

cd pre-solver
./clean-cmake.sh
cd -
