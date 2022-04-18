#!/bin/bash

valgrind --tool=cachegrind --I1=32768,2,64 --D1=32768,2,64 --L2=4194304,8,64 $@

