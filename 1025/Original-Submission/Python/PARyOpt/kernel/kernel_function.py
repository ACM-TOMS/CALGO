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

## --- end license text --- ##
import numpy as np


class KernelFunction:
    """
    Attributes:

    :var theta: Kernel parameters. Can be None, a scalar, or a vector. Should be initialized to the appropriate \
    type/len during initialization.
    :var theta0: Scaling factor
    """

    theta = None
    theta0 = 1.0

    def __init__(self):
        pass

    def eval(self, x1: np.array, x2: np.array) -> float:
        """
        Evaluate the kernel function

        :param x1: location 1
        :param x2: location 2
        :return:
        """
        raise NotImplemented('Kernel function eval not implemented')

    def derivative(self, x1: np.array, x2: np.array) -> np.array:
        """
        Evaluate the derivative of the kernel

        :param x1: location 1
        :param x2: location 2
        :return:
        """
        raise NotImplemented('Kernel function derivative not implemented')

