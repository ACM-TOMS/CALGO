 
   Table of Contents
   -----------------

   1. Installation

   2. Description

   3. Execution of MULTIPLICITY

   4. Files and the sub-folder 

   5. Disclaimer

   6. Contact information


1. Installation
   ------------

   Step 1. Unzip the file "Multiplicity.zip", generating the folder 
               "Multiplicity"

   Step 2. Set the Matlab path to the folder "Multiplicity"
             e.g.   >> path(path,'c:/Test/Multiplicity')

2. Description
   ------------

The m-file "multiplicity.m" computes the mulitiplicity structure of a
nonlinear system at a numerical zero. For details, see 
      An algorithm and software for computing multiplicity structures at
      zeros of nonlinear systems, W. Hao, A. J. Sommesse and Z. Zeng 


3. Execution of MULTIPLICITY
   ------------

To compute the multiplicity, the user needs to set up input items including
the nonlinear system, the variable names, and the zero. 

For example, consider the nonlinear system

               sin(x)*cos(y)-x       = 0
               sin(y)*sin(x)^2 - y^2 = 0

at the zero (0,0), the multiplicity can be computed by the following
statements:

    >> f = {'sin(x)*cos(y)-x', 'sin(y)*sin(x)^2-y^2'}; % set the system
    >> variables = {'x','y'};                          % set the variables
    >> zero = [0,0];                                   % set the zero
    >> m = multiplicity(f,variables,zero)         % compute the multiplicity

The full-featured call can be carried out by

  >> [m, D, H] = multiplicity(f, variables, zero, options)

where the output D contains a basis for the dual space, and H is the Hilbert
function value. See MULTIPLICITY for details.

Use "help multiplicity" for the inline documentation. 

4. Files and the sub-folder 
   ------------

   multiplicity.m:  The main code

   Demo.m:          Matlab script for a demostration

   Readme.txt:      This file
 
   TestSuite:       The folder containing the benchmark nonlinear systems 
                    for testing. See the file"TestReadme.txt" in the folder.

5. Disclaimer
   ------------

MULTIPLICITY is distributed free of charge on an "as is" basis. Its 
intended usage is educational, so that the user may gain a greater 
understanding of multiplicity structure of nonlinear systems. Any other use 
is strictly the user's responsibility.  In publications based on results 
obtained using MULTIPLICITY or its successors, the use of MULTIPLICITY 
should be acknowledged. Any distribution of derived codes that extend or 
modify MULTIPLICITY should acknowledge the original source and authorship.

6. Contact Information
   ------------

If you publish a result using MULTIPLICITY, we would like to hear about it. 
You may write to one of us, Wenrui Hao (whao@nd.edu), Andrew Sommese
(sommese@nd.edu) and Zhonggang Zeng(z-zeng@neiu.edu).

------------------------------------------------------------------------
 MULTIPLICITY Version 1.0 Date: 2011/12/02 
------------------------------------------------------------------------