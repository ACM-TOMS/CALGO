from distutils.core import setup, Extension
from Cython.Build import cythonize

import os
import sys


    
if (len(sys.argv) != 4 and len(sys.argv) != 2):
    print("The script ", str(sys.argv[0]), " requires 3 arguments (building in place) or 1 arguments (installing). When building in place, arguments 1 is for passing Nomad path.")
    exit()

root_build_dir=""
os_include_dirs=""
if (len(sys.argv) == 4):
    #print("original sys argv: " + str(sys.argv))
    root_build_dir = str(sys.argv[1])   # Argument 1 is to path to find libraries and headers
    os_include_dirs = [root_build_dir + "/../../src"]
    del sys.argv[1]
    # print("new sys argv: " + str(sys.argv))

build_lib_dir = root_build_dir + "/src"
installed_lib_dir1 = root_build_dir + "/lib"
installed_lib_dir2 = root_build_dir + "/lib64"

compile_args = []
compile_args.append("-w")
#compile_args.append("-Wall")
if not sys.platform.startswith("win"):
    compile_args.append("-std=c++14")
    # compile_args.append("-Wextra")
    compile_args.append("-pthread")

link_args = []

# Use gcc
if (str(os.environ.get("CC")) == ""):
   os.environ["CC"] = "gcc"
if (str(os.environ.get("CXX")) == ""):
   os.environ["CXX"] = "g++"

# Look for librairies in Nomad distribution
# Default variables are for Linux
lib_extension="so"

# The rpath is set to the libraries installation directory.
if not sys.platform.startswith("win"):
    link_args.append("-Wl,-rpath," + installed_lib_dir1)
    link_args.append("-Wl,-rpath," + installed_lib_dir2)

if sys.platform == "darwin":
     link_args.append("-headerpad_max_install_names")
     lib_extension = "dylib"
elif sys.platform.startswith("win"):
     lib_extension = "lib"

# For building/linking the interface, the libraries build directory is used.
if sys.platform.startswith("win"):
    link_args.append(build_lib_dir + "/Release/nomadUtils." + lib_extension)
    link_args.append(build_lib_dir + "/Release/nomadEval." + lib_extension)
    link_args.append(build_lib_dir+ "/Release/nomadAlgos." + lib_extension)
else:
    link_args.append(build_lib_dir + "/libnomadUtils." + lib_extension)
    link_args.append(build_lib_dir + "/libnomadEval." + lib_extension)
    link_args.append(build_lib_dir+ "/libnomadAlgos." + lib_extension)

setup(
    name='PyNomad',
    ext_modules = cythonize(Extension(
           "PyNomad", # extension name
           sources = ["PyNomad.pyx"], # Cython source and interface
           include_dirs = os_include_dirs,
           extra_compile_args = compile_args,
           extra_link_args = link_args,
           language = "c++"
           ))
)

