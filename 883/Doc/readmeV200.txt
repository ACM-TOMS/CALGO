
Revision of SparsePOP200 from SparsePOP120
1) C++ subroutines are added to speedup the construction of asparse SDP relaxation problem. Two kinds of SparsePOP are now included:A set of MATLAB programs combined with C++ subroutines and a set of all 
MATLAB programs. Input and output are identical for the both SparsePOPs;they construct the same SDP relaxation problem.2) User interface of SparsePOP is improved so that it can now accept
three types of inputs: the name of a file describing a POP in the GAMS
scalar format, the name of a file describing a POP in the SparsePOP format,
and a set of input arguments for a POP.3) Some bugs are fixed.