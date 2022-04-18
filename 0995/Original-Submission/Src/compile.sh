#!/bin/bash

#The build script will use the compiler that gcc and g++ point to

cd src/triangle
make triangle.o
cd ../
cmake . .
make