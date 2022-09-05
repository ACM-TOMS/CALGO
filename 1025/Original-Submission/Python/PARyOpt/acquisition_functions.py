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
"""
Contains a library of acquisition functions that can be used in bayesian optimization
"""
import numpy as np

import PARyOpt.utils as utils


def lower_confidence_bound(mean: float, variance: float, curr_best: float = 0., kappa: float = 1.) -> float:
    """
    lower confidence bound of improvement : used for minimization problems

    :param mean: mean of surrogate
    :param variance: variance of surrogate
    :param curr_best: current best evaluated point
    :param kappa: exploration - exploitation tradeoff parameter
    :return: lower confidence bound
    """
    return mean - kappa * variance - curr_best


def upper_confidence_bound(mean: float, variance: float, curr_best: float=0., kappa: float = 1.) -> float:
    """
    upper confidence bound of improvement: used in the case of maximization problems

    :param mean: mean of surrogate
    :param variance: variance of surrogate
    :param curr_best: current best evaluated point
    :param kappa: exploration - exploitation tradeoff parameter
    :return: upper confidence bound
    """
    return mean + kappa * variance - curr_best


def expected_improvement(mean: float, variance: float, curr_best: float=0., _: float = 1.) -> float:
    """
    Expected improvement of objective function \
    'A Tutorial on Bayesian Optimization of Expensive Cost Functions, \
    with Application to Active User Modeling and Hierarchical Reinforcement Learning'

    :param mean: mean of surrogate
    :param variance: variance of surrogate
    :param curr_best: current best evaluated point
    :param kappa: exploration - exploitation tradeoff parameter
    :return: expectation of improvement
    """
    if variance > 1e-30:
        gamma = -1.0 * (curr_best - mean) / np.sqrt(variance)
        return np.sqrt(variance) * (gamma * utils.cdf_normal(gamma) + utils.pdf_normal(gamma))
    else:
        return 0.0


def probability_improvement(mean: float, variance: float, curr_best: float=0., _: float = 1.) -> float:
    """
    Probability of improvement of objective function \
    'A Tutorial on Bayesian Optimization of Expensive Cost Functions, \
    with Application to Active User Modeling and Hierarchical Reinforcement Learning'

    :param mean: mean of surrogate
    :param variance: variance of surrogate
    :param curr_best: current best evaluated point
    :param kappa: exploration - exploitation tradeoff parameter
    :return: probability of improvement
    """
    if variance > 1e-16:
        gamma = -(curr_best - mean) / np.sqrt(variance)
        return utils.cdf_normal(gamma)
    else:
        return 0.0
