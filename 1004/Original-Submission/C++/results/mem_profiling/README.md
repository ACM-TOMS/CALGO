
Memory profiling
----------------

The standalone program defined in mem_profile.cpp calculates a single log signature.
It is built with `mem_profile_build.sh`. `mem_profile_run.py` uses the `massif` tool
from `valgrind` to work out the total memory usage. `mem_profile_analyse.py` plots
results.