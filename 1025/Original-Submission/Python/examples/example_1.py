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
"""
Example 1: Using user-defined functions for initialization, kernel, acquisition and acquisition optimization
1- dimensional
We shall solve the same problem as example 0 with the same internal functions but supplied by the user

What you will learn:
* user-defined initialization
* user-defined kernel function
* user-defined acquisition function
* user-defined acquisition optimizer
"""
import logging
import time

import numpy as np
from scipy import optimize

import examples.examples_all_functions as exf
from typing import List, Callable
from PARyOpt import BayesOpt, utils
from PARyOpt.evaluators import FunctionEvaluator

from PARyOpt.kernel import KernelFunction


# custom initialization
def my_init_strategy(n_dim: int, n_init: int, l_bound: np.array, u_bound: np.array) -> List[np.array]:
    """
    User customized initialization -- since the example is one-dimensional
    we shall uniformly divide the domain into n_init parts and then return the values
    :param n_dim: dimensionality
    :param n_init: number of initial points
    :param l_bound: lower bound
    :param u_bound: upper bound
    :return:
    """
    if n_dim == 1:
        lin_spaced_intervals = np.linspace(l_bound, u_bound, n_init+2)
        return [np.asarray([lin_spaced_intervals[i+1]]) for i in range(n_init)]
    else:
        raise NotImplementedError("init strategy not implemented for more than 1 dimension")


# custom kernel function
class MyKernelFunction(KernelFunction):
    """
    user customized kernel function
    """

    theta = 1.0

    def eval(self, x1: np.array, x2: np.array) -> float:
        """
        actual kernel function
        :param x1:
        :param x2:
        :return:
        """
        x1 = np.divide(x1, self.theta)
        x2 = np.divide(x2, self.theta)
        dist = utils.distance(x1, x2)
        rval = np.sqrt(3.0) * dist
        return self.theta0 * (1+rval) * np.exp(-rval)

    def derivative(self, x1: np.array, x2: np.array) -> np.array:
        """
        derivative of kernel function
        currently not useful, so we will not implement anything here
        :param x1:
        :param x2:
        :return:
        """
        pass


# custom acquisition function
def my_acquisition_function(mean: float, variance: float, curr_best: float = 0., kappa: float = 1.) -> float:
    """
    user customized acquisition function -- negative log of LCB
    :param mean:
    :param variance:
    :param curr_best:
    :param kappa:
    :return:
    """
    return -1.0 * float(np.log(abs(mean - kappa * variance - curr_best)))


# custom acquisition function optimizer
def my_acq_optimizer(func: Callable[[np.array], float], x0: np.array, l_bound: np.array, u_bound: np.array) -> np.array:
    """
    user customized acquisition function optimizer
    same as the default optimizer, i.e., POWELL search from scipy.optimize
    :param func:
    :param x0:
    :param l_bound:
    :param u_bound:
    :return:
    """
    bounds = list(zip(l_bound, u_bound))
    res = optimize.minimize(fun=func, x0=x0, bounds=bounds, method='Powell')
    return res.x


if __name__ == "__main__":
    # logging setup -- if this is not done, it will be streamed to stdout
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)  # either NOTSET, INFO, DEBUG, WARNING, ERROR, CRITICAL -- different levels of log
    log_file_name = '../logs/example1_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S"))
    fh = logging.FileHandler(log_file_name, mode='a')

    # log file format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # first line of the log
    logger.info('Example run for Bayesian optimization')

    # let us define basic parameters here

    # dimensionality of problem
    n_dim = 1

    # bounds, as a numpy array
    # l_bound[0] is the lower bound for dimension 0, etc.
    l_bound = np.asarray([-12.])
    u_bound = np.asarray([12.])

    # number of cost function calls to create the initial surrogate
    n_init = 2

    # maximum number of iterations
    max_iter = 10

    # cost function : simple parabola
    cost_function = exf.parabolic_cost_function
    # how to evaluate the cost function - we will use the serial, in-script evaluator
    evaluator = FunctionEvaluator(cost_function)

    # kappa strategy : defines exploration vs exploitation .. we will use a pure exploration here
    # large values of kappa are exploratory (prioritize areas where confidence is low),
    # small values are exploitatory (prioritize areas where the surrogate seems optimal)
    # the following definition returns a value of 1000 for any iteration value
    def my_kappa(iteration: int) -> float:
        """
        return a constant kappa value of 1000.
        :param iteration: iteration
        :return:
        """
        return 1000.0


    # instantiate the optimizer with cost function, bounds, initial evaluation points,
    # type of kernel function (squared exponential), acquisition function (lower confidence bound)
    # and a kappa strategy (constant here)

    b_opt = BayesOpt(cost_function=evaluator,
                     l_bound=l_bound, u_bound=u_bound, n_dim=n_dim,
                     n_init=2, init_strategy=my_init_strategy,
                     kern_function=MyKernelFunction(),
                     acq_func=my_acquisition_function,
                     kappa_strategy=my_kappa,
                     acq_func_optimizer=my_acq_optimizer,
                     if_restart=False)

    logger.info('BO initialized')

    # visualize the surrogate generated by the optimizer, with no hyper-parameter optimization
    exf.visualize_fit(b_opt)
    # estimate best kernel parameters: perform hyper parameter optimization
    b_opt.estimate_best_kernel_parameters(theta_bounds=[[0.1, 10]])
    exf.visualize_fit(b_opt)

    # an example of how to evaluate the surrogate at any point and see how far it is from the actual cost function
    pt_to_evaluate = np.asarray([0.0])
    mean, variance = b_opt.evaluate_surrogate_at(pt_to_evaluate)
    true_cost_function = exf.parabolic_cost_function(pt_to_evaluate)

    print('At {}, surrogate mean:{}, variance:{}, true value:{}'.format(pt_to_evaluate, mean, variance,
                                                                        true_cost_function))

    # update iterations
    for curr_iter in range(max_iter):
        # update_iter finds acquisition function optima and updates the prior to get the posterior
        # (optimize acquisition function to pick the next set of points,
        # evaluate the cost function at those points, and update the surrogate with the new values)
        b_opt.update_iter()
        # estimate and set best kernel parameters
        b_opt.estimate_best_kernel_parameters(theta_bounds=[[0.1, 10]])
        # visualize it every other iteration
        if curr_iter % 2:
            exf.visualize_fit(b_opt)

    # get the population and respective function values from the optimizer
    total_population, function_values = b_opt.get_total_population()
    # get current best evaluated value
    best_location, best_value = b_opt.get_current_best()

    result_txt = 'Optimization done for {} iterations, best evaluation is at {} with cost: {}'. \
        format(b_opt.get_current_iteration(), best_location, best_value)

    logger.info(result_txt)
    print(result_txt)

    csv_file = '../temp/example_1_data.csv'
    b_opt.export_csv(csv_file)
    logger.info('Data exported to {}'.format(csv_file))
