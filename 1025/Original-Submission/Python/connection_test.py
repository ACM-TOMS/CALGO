"""
    Copyright (c) 2018 Baskar Ganapathysubramanian, Balaji Sesha Sarath Pokuri
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
"""

## --- end license text --- ##
from PARyOpt.evaluators.connection import Connection, Host
from PARyOpt.evaluators.async_sbatch import VALUE_FROM_FILE, AsyncSbatchEvaluator

import numpy as np
import os
import stat
import getpass


def job_gen(lcl_dir: str, x: np.array) -> None:
    path = os.path.join(lcl_dir, 'test.py')
    with open(path, 'w') as f:
        f.write('#!/bin/env python\n')
        f.write('print(123456.789)\n')
    os.chmod(path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)


job_script = """#!/bin/bash
#SBATCH --job-name='paryopt_test'
#SBATCH --output='output.txt'
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:02:00

module load python
python test.py > value.txt
"""

host = Host(getpass.getuser(), 'condo.its.iastate.edu')
# host = Host(getpass.getuser(), 'stampede2.tacc.utexas.edu')
funcEval = AsyncSbatchEvaluator(host, job_gen, job_script, lcl_parse_result=VALUE_FROM_FILE('value.txt'))

x = np.array([1.0])
xs = [x]
completed, pending, failed = funcEval.evaluate_population(xs)
print([str(c) for c in completed])
print([str(f) for f in failed])