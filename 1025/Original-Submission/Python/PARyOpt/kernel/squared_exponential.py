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
from .kernel_function import KernelFunction
from ..utils import distance


class SquaredExponential(KernelFunction):
    """
    Squared exponential rbf
    """

    theta = 1.0

    def eval(self, x1: np.array, x2: np.array) -> float:
        assert (not hasattr(self.theta, '__len__')) or x1.shape == self.theta.shape
        x1 = np.divide(x1, self.theta)
        x2 = np.divide(x2, self.theta)
        dist = distance(x1, x2)
        return self.theta0 * np.exp(-0.5 * dist * dist)

    def derivative(self, x1: np.array, x2: np.array) -> np.array:
        assert (not hasattr(self.theta, '__len__')) or x1.shape == self.theta.shape
        x1 = np.divide(x1, self.theta)
        x2 = np.divide(x2, self.theta)
        dist = distance(x1, x2)
        return -self.theta0 * np.exp(-0.5 * dist * dist) * (x1 - x2) / self.theta
