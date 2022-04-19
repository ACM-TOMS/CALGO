
This example illustrates how to use NOMAD with a MATLAB black-box problem.

MATLAB must have been installed with the MCC compiler, in order to create a stand-alone
executable from the example program file bb.m.
From the MATLAB command line, type 'mcc -m bb' to generate the executable file bb.exe.

From the command prompt, you can test the validity of the black-box with the command 'bb.exe x.txt'.

One the black-box is tested, use nomad with the parameters file parameters.txt.