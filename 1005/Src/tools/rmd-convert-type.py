#!/usr/bin/python
""" A simple program to convert BLAS-RMD programs from single precision
to double precision.

Usage: rmd-convert-type.py [input_file] [-o output_file]

Options and arguments:
input_file     : Name of file to be converted. If `no input_file` is given,
                 then all input will be taken from the standard input.
-o output_file : Write output into output_file Warning: Any file with
                 name `output_file` will be overwritten. If no `output_file` 
                 is given, then all output will be written to the standard
                 output.
"""

import sys
import getopt
import re

## Various Fortran numerical values used when converting
## BLAS operations to be converted
# generic names
BLAS_OPS = [
    'swap', 'scal', 'copy', 'axpy', 'dot', \
    'gemv', 'gbmv', 'symv', 'sbmv', 'spmv', \
    'trmv', 'tbmv', 'tpmv', 'trsv', 'tbsv', \
    'tpsv', 'ger', 'syr', 'spr', 'syr2', \
    'spr2', 'gemm', 'symm', 'syrk', 'syr2k', \
    'trmm', 'trsm', 'potrf', 'nrm2', 'rot', \
    'rotg', 'rotm', 'rotmg', 'asum']
BLAS_OPS  = list(map(lambda x: 's'+x, BLAS_OPS))
BLAS_OPSC = list(map(lambda x: x.upper(), BLAS_OPS))

## Conversion functions
# A conversion function from prefix `x` to `y` should replace every `xOP`
# from `xBLAS_OPS` found in the input string with `yOP` and use
# regex matching to replace `rmd_x*` with the corresponding `rmd_y*`.
# It should also replace type declarations accordingly, for example
# `real` becomes `double precision`.

def num_convert(matchobj):
    result = matchobj.group(1) + matchobj.group(2) + 'd'
    if matchobj.group(4):
        result += matchobj.group(4)
    else:
        result += '0'
    return result

def s2d(input_string):
    # Mutable list for in-place replacing of certain values
    output_builder = list(input_string)
    for op in BLAS_OPS: # Change prefix of blas and BLAS-RMD operations
        for match in re.finditer(op, input_string):
            output_builder[match.start()] = 'd'
    for op in BLAS_OPSC: # For upper case blas_ops
        for match in re.finditer(op, input_string):
            output_builder[match.start()] = 'D'
    # Change prefix for utility functions
    output_string = str().join(output_builder) # Prepare output string
    output_string = output_string.replace('sdddot', 'sdsdot') # Change back
    output_string = output_string.replace('rmd_s', 'rmd_d') # Prefix for utility fns
    output_string = output_string.replace('RMD_S', 'RMD_D') # Prefix for utility fns
    output_string = re.sub(r'\breal\b', 'double precision', output_string) # Replace type
    output_string = re.sub(r'\bREAL\b', 'double precision', output_string) # Replace type
    output_string = re.sub(r'\bsngl\b', 'dble', output_string) # Replace type
    num_pattern = r'([^a-zA-Z0-9])(\d+\.\d*)([eE]([+-]?\d+))?' # Replace literals
    num_pattern1 = r'([^a-zA-Z0-9])(\d+)([eE]([+-]?\d+))' # Replace literals
    # 5.32 --> 5.32d0, 5.32e9 --> 5.32d9
    output_string = re.sub(num_pattern, num_convert, output_string)
    # 5e9 --> 5d9
    output_string = re.sub(num_pattern1, num_convert, output_string)
    output_string = output_string.replace('isdoubleversion = .false.',
                                          'isdoubleversion = .true.')
    return output_string

def main(argv=None):
    # Note: Currently, the whole stream/file is stored as a string in
    #       memory. For very long streams/files, it might be
    #       better to iterate on a per-line basis.
    
    # Get command line arguments
    if argv == None:
        argv = sys.argv[1:]
    options, remainder = getopt.gnu_getopt(argv, "o:")
    
    # Read input
    f=open(remainder[0], 'r')
    input_string = f.read()
    f.close()
    
    # Process input
    output_string = s2d(input_string)
    
    # Write output
    if options == []:
        sys.stdout.write(output_string)
    else:
        for opt, arg in options:
            if opt == '-o':
                f=open(arg, 'w')
                f.write(output_string)
                f.close()
    
if __name__ == "__main__":
    main()
