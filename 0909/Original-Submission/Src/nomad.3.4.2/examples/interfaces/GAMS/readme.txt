
This example illustrates how to use NOMAD with a GAMS black-box problem.

The example runs on windows.

The GAMS program is located in file bb.gms.

The source code is in bb.cpp: it is a C++ wrapper that calls GAMS.
It can be compiled with visual C++ or minGW.

The wrapper can be tested with the command 'bb.exe .\points\x0.txt'.

Once the wrapper is compiled, use nomad with the parameters file parameters.txt.
