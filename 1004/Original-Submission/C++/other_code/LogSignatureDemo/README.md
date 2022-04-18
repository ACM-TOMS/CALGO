##Demonstration of Log Signature calculations

This is an example of calculating the log signature of a path using the method described [here](http://www2.warwick.ac.uk/fac/cross_fac/complexity/people/students/dtc/students2013/reizenstein/logsignatures.pdf).

The file `logsignature.py` is the same as the original distribution on [Jeremy's website](http://www2.warwick.ac.uk/fac/cross_fac/complexity/people/students/dtc/students2013/reizenstein) except you can use command-line arguments to specify the dimensions (D) and levels (M). Each argument is a comma separated list of ranges or values.

###Usage
First save [this file from Casas and Murua](http://www.ehu.eus/ccwmuura/research/bchLyndon20.dat) to this directory.
To output the log signature of a random path in 3D up to level 6:

```
python logsignature.py 3 6
g++ --std=c++11 -DDIM=3 -DLEVEL=6 bch.cpp calculate.cpp
./a.out 
```

To generate the BCH code for 2 and 3 dimensions for levels 2 through 7:

```
python logsignature.py 2,3 2-7
```

Note that running the python script (which generates the files `bch.h` and `bch.cpp`) and compiling the code will in general be slow, but the final program should be fast. The generated code file will in general be large, as will its associated object file if you generate it. The compiling can be quicker with clang++ than g++. 

Works with: python 2.6+ or 3.x .
