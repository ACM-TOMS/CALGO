"""
---
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
---
"""
import sys
import numpy as np
import time
import examples.examples_all_functions as exf

# read the command line arguments into an array
xs = np.asarray([float(x) for x in sys.argv[1:]])

# evaluate the cost, i.e. the parabola
cost = exf.parabolic_cost_function(x=xs)

# In order to demonstrate uncertain evaluation times, we shall use a random sleep in each cost function.
# In this case, each function evaluation can take between 0-10 sec
time.sleep(np.random.random() * 10.)

# Write into a result file. Note that this script is evaluated in its respective folder and so
# the result.txt file will be in the generated folder and not the home directory of running example3.py
with open('result.txt', 'w') as f:
    f.write(str(cost) + '\n')
