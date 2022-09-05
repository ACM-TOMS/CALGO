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
Collection of functions for the examples: plug in different functions to get each example
"""

import numpy as np
import scipy
import scipy.io as io

import logging
import time
import matplotlib.pyplot as plt
import csv
from typing import Callable, List
from matplotlib import cm, gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize

from PARyOpt import BayesOpt
from PARyOpt.evaluators import FunctionEvaluator
import PARyOpt.utils as ut


def plot_only_vertical_lines_bayes_opt_population(b_opt: BayesOpt) -> None:
    """
    Plots vertical lines at the locations of function evaluation
    Creates a plot with only the function evaluation
    :param b_opt:
    :return:
    """
    x_eval, y_eval = b_opt.get_total_population()
    y_eval = np.asarray(y_eval)
    plt.plot(x_eval, y_eval, 'r*')
    for x_now, y_now in zip(x_eval, y_eval):
        plt.plot([x_now, x_now], [y_now, 0.0], 'b--', linewidth=2)
    plt.plot([b_opt.l_bound, b_opt.u_bound], [0.0, 0.0], '--', color='grey')
    plt.title('Number of evaluation points: {}'.format(b_opt.curr_points))
    plt.show()


def visualize_fit(b_opt: BayesOpt) -> None:
    """
    Visualize the 2d fit
    :param b_opt: BayesOpt object
    :return: None
    """
    if b_opt.n_dim == 1:
        xs = np.linspace(start=b_opt.l_bound, stop=b_opt.u_bound, num=1000, retstep=False)
        _mean = np.empty(xs.shape[0])
        _variance = np.empty(xs.shape[0])
        _acq = np.empty(xs.shape[0])
        _penalty = np.empty(xs.shape[0])
        _s_opt = np.empty(xs.shape[0])
        kappa = b_opt.get_kappa[0](b_opt.curr_iter+1)
        for _i in range(xs.shape[0]):
            # _target[_i] = cost_function(xs[_i])
            _mean[_i], var = b_opt.evaluate_surrogate_at(np.asarray([xs[_i]]))
            _variance[_i] = np.sqrt(abs(var))
            _acq[_i] = b_opt.acquisition_function[0](_mean[_i], _variance[_i])
            _penalty[_i] = b_opt.penalty_function(np.asarray([xs[_i]]), kappa)
            _s_opt[_i] = b_opt.acquisition(np.asarray([xs[_i]]), kappa)
            min_s = np.argmin(_s_opt)
        # gather the evaluation points
        x_eval, y_eval = b_opt.get_total_population()
        y_eval = np.asarray(y_eval)
        gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
        plt.title('Number of evaluation points: {}'.format(b_opt.curr_points))
        ax0 = plt.subplot(gs[0])
        ax0.plot(x_eval, y_eval, 'r*', linewidth=4, label='Actual')
        ax0.plot(xs, _mean, 'g-', linewidth=2, label='Surrogate mean')
        ax0.plot([b_opt.l_bound, b_opt.u_bound], [0.0, 0.0], '--', color='grey')
        for x_now, y_now in zip(x_eval, y_eval):
            ax0.plot([x_now, x_now], [y_now, 0.0], 'y--', linewidth=3)
        ax1 = plt.subplot(gs[1])
        ax1.plot(xs, _variance, '-', linewidth=3, label='Variance')
        ax1.plot([b_opt.l_bound, b_opt.u_bound], [0.0, 0.0], '--', color='grey')
        plt.fill_between(xs, _mean+_variance, _mean-_variance, alpha=0.15)
        plt.plot(xs, _mean + _variance, 'b--', label='mean+variance', alpha=0.2)
        plt.plot(xs, _mean - _variance, 'k--', label='mean-variance', alpha=0.2)
        plt.legend(loc='upper left')
        plt.show()

    elif b_opt.n_dim == 2:
        print('2D visualization!')
        xs = np.linspace(start=b_opt.l_bound[0], stop=b_opt.u_bound[0], num=50, retstep=False)
        ys = np.linspace(start=b_opt.l_bound[1], stop=b_opt.u_bound[1], num=50, retstep=False)
        _mean = np.empty((xs.shape[0], ys.shape[0]))
        _variance = np.empty((xs.shape[0], ys.shape[0]))
        X, Y = np.meshgrid(xs, ys)
        for _i, _j in [(x, y) for x in range(xs.shape[0]) for y in range(ys.shape[0])]:
            location = np.asarray([xs[_i], ys[_j]])
            _mean[_i, _j], var = b_opt.evaluate_surrogate_at(location)
            _variance[_i, _j] = np.sqrt(abs(var))
        tot_pop, func_pop = b_opt.get_total_population()

        xs = [p[0] for p in tot_pop]
        ys = [p[1] for p in tot_pop]

        vmax = max(np.array(func_pop))
        vmin = min(np.array(func_pop))
        plt.figure(1)
        ax1 = plt.subplot(1, 2, 1)
        plt.plot(xs, ys, 'r*')
        ax1.set_aspect('equal')
        ax1.set_title('Mean of estimate')
        p1 = plt.contourf(X, Y, _mean.transpose(), cmap=cm.binary, vmin=vmin, vmax=vmax)
        plt.plot(xs, ys, 'r*')

        ax2 = plt.subplot(1, 2, 2)
        ax2.set_aspect('equal')
        ax2.set_title('Variance of estimate')
        p2 = plt.contourf(X, Y, _variance.transpose(), cmap=cm.binary)
        plt.plot(xs, ys, 'r*')
        plt.show(block=True)
        # plot wireframe plots
        fig2 = plt.figure(2)
        ax1 = fig2.add_subplot(1, 2, 1, projection='3d')
        ax1.plot_wireframe(X, Y, _mean.transpose())
        ax1.scatter(xs, ys, func_pop, s=14, c='red')
        ax1.set_aspect('equal')
        ax1.set_title('Mean of estimate')
        ax2 = fig2.add_subplot(1, 2, 2, projection='3d')
        ax2.plot_wireframe(X, Y, _variance.transpose())
        ax2.set_aspect('equal')
        ax2.set_title('Variance of estimate')
        plt.show(block=True)

    else:
        raise NotImplemented()


def add_points_from_file(bo: BayesOpt, filename: str) -> BayesOpt:
    """
    Adds points from the file filename and returns back the bayes opt object
    1D kriging
    :param filename: file that contains the data in a row format
    :param bo: bayes opt object
    :return: updated bayes opt object
    """
    with open(filename, mode='r') as csvfile:
        csv_file_lines = csv.reader(csvfile, delimiter=',')
        for row_num, row in enumerate(csv_file_lines):
            if row_num == 0:
                # skipping the header
                pass
            else:
                bo.add_point(x=np.asarray([float(row[0])]), y=float(row[1]),
                             if_check_nearness=True)
    bo.update_surrogate()
    return bo


def custom_sum_constraint(phi: np.array) -> bool:
    """
    :param phi: 
    :return: 
    """
    return phi.sum() <= 100.0


def my_kappa_annealed(curr_iter: int, t_const: float = 0.5, amp: float = 50) -> float:
    """
    Annealed kappa
    :param amp: starting amplitude
    :param curr_iter: current iteration
    :param t_const: time constant for exponential
    :return: value of kappa
    """
    kappa = amp * np.exp(-t_const*curr_iter)
    return float(kappa)


def my_kappa_sine(curr_iter: int, t_period: float = 2*np.pi, amp: float = 1.0 ) -> float:
    """
    Sine variation of kappa
    :param curr_iter: current iteration
    :param t_period: time period of sine curve
    :param amp : amplitude
    :return:
    """
    kappa = amp * np.sin(curr_iter * t_period)
    return float(kappa)


def my_kappa_constant(curr_iter: int, kappa_const: float = 1.0)-> float:
    """
    returns a constant kappa value
    :param kappa_const: value to return
    :param curr_iter: current iteration
    :return:
    """
    return float(kappa_const)


def ackley_function(x: np.array) -> float:
    """
    Ackley function
    test area : -32.768 <= x[i] <= 32.768
    minima at x[i] = 0.  : f(x) = 0.0
    :param x:
    :return:
    """
    sum1 = np.mean(x**2)
    sum2 = np.mean(np.cos(2*np.pi*x))
    return float(-20.0 * np.exp(-0.2 * np.sqrt(sum1)) - np.exp(sum2) + 20 + np.exp(1.0))


def rosenbrock(x: np.array) -> float:
    """
    Rosenbrock function
    test area : -2.048 <= x[i] <= 2.048
    minima at x[i] = 1.0  : f(x) = 0.0
    :param x:
    :return:
    """
    rosenbrock_sum = 0
    dim = x.shape[0]
    for i in range(dim - 1):
        rosenbrock_sum += (100 * (x[i+1] - x[i]**2)**2 + (1.0 - x[i])**2)

    return float(rosenbrock_sum)


def rastrigin(x: np.array) -> float:
    """
    Rastrigin function
    test area : -5.12<=x[i]<=5.12
    minima at x[i] = 0.0 : f(x) = 0.0
    :param x:
    :return:
    """
    dim = x.shape[0]
    rastrigin_sum = 0.0
    for i in range(dim):
        rastrigin_sum += (x[i]**2 - 10.0 * np.cos(2*np.pi*x[i]))
    return float(10*dim + rastrigin_sum)


def griewangk(x: np.array) -> float:
    """
    Griewangk's function
    test area: -600 <= x[i] <= 600
    minima at x[i] = 0.0 : f(x) = 0.0
    :param x:
    :return:
    """
    prod = 1.0
    dim = x.shape[0]
    for i in range(dim):
        prod *= np.cos(x[i]/np.sqrt(i+1))
    return float(1./4000 * np.sum(x**2) - prod + 1.0)


def constant_cost_function(x: np.asarray) -> float:
    """
    does nothing
    :param x:
    :return: constant value of 100.0
    """
    return 100.0


def parabolic_cost_function(x: np.array) -> float:
    """
    y = (x-2.5) ^ 2 + 5
    Dimension independant
    :param x: location
    :return:
    """
    y = np.sum((x-2.5) ** 2 + 5)
    return float(y)


def multi_start_powell_acq_optimizer(func: Callable[[np.array], float], x0: np.array, l_bound: np.array, u_bound:np.array) -> np.array:
    """
    Multi start local optimizer, using powell method
    uses latin hypercube sampling for starting points
    :param func: acquisition function to be minimized
    :param x0: initial guess
    :param l_bound: lower bound
    :param u_bound: upper bound
    :return: minima of func
    """
    n_dim = x0.shape[0] # type: int
    starting_point_list = ut.lhs(n_dim, samples=6 * n_dim, criterion='cm', iterations=20)   # type: List([np.array])
    optima_list = []    # type: list[np.array]
    optima_value_list = []  # type: list[float]
    # find the respective minima from the various starting points
    for pt in starting_point_list:
        bounds = list(zip(l_bound, u_bound))
        res = optimize.minimize(fun=func, x0=x0, bounds=bounds, method='Powell')
        optima_list.append(res.x)
        optima_value_list.append(res.func)
    # find the smallest of the found minima
    min = optima_value_list[0]
    min_point = optima_list[0]
    for opt, opt_val in zip(optima_list, optima_value_list):
        if opt_val < min:
            min = opt_val
            min_point = opt
    # return the (global) minima from the available local optima
    return min_point


def my_acq_optimizer_brute(func: Callable[[np.array], float], x0: np.array, l_bound: np.array, u_bound:np.array) -> np.array:
    """
    brute search algorithm for acquisition optimization
    :param u_bound: upper bound
    :param l_bound: lower bound
    :param func: funciton to minimize
    :param x0: Initial guess
    :return:
    """
    bounds = list(zip(l_bound, u_bound))
    res = scipy.optimize.brute(func, bounds, Ns=50, full_output=True, finish=None)
    return res[0]


if __name__ == "__main__":

    # logging setup -- if this is not done, it will be streamed to stdout
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)   # can be NOTSET, INFO, DEBUG, WARN, ERROR, CRITICAL
    # log file name
    fh = logging.FileHandler('examples_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S")), mode='a')

    # log file format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # first line of the log
    logger.info('Example run for Bayesian optimization')

    # Local, in-script evaluation
    # this can be either 
    # >>> FunctionEvaluator 
    # >>> AsyncLocalEvaluator
    # >>> AsyncSbatchEvaluator
    # more examples with different evaluators will be added soon
    evaluator = FunctionEvaluator(constant_cost_function)

    # set up problem related parameters
    n_dim = 2                               # dimensionality
    l_bound = np.asarray([0.0, 0.0])        # lower bound
    u_bound = np.asarray([100.0, 100.0])    # upper bound
    n_opt = 1                               # max number of optima per iteration
    n_init = 2                              # initial sampling for Latin hypercube sampling

    # Initialize bayesian optimization
    # takes in cost function, dimensionality, bounds, optima per iteration, initialization parameters
    # kernel function information, acquisition function info, kappa strategy, restart info and custom constraints
    bo = BayesOpt(cost_function=evaluator, n_dim=n_dim, n_opt=n_opt,
                  n_init=n_init, kern_function='matern_52', init_strategy=None, u_bound=u_bound,
                  l_bound=l_bound, acq_func='LCB', kappa_strategy=lambda _: 10000.0, if_restart=False,
                  do_init=True, constraints=[custom_sum_constraint])
    # This step sets up the object and does `n_init` function evaluations serially on the local machine 

    logger.info('Bayes Opt initialized')

    # user can add points manually from a file:
    manual_data_to_enter_file= '/home/balajip/Desktop/bo_md/ternary_mixing.csv'
    bo = add_points_from_file(bo, filename=manual_data_to_enter_file)
    # update the surrogate
    bo.update_surrogate()
    # visualize
    visualize_fit(bo)

    logger.info('Bayes Opt -- data added from file: {}'.format(manual_data_to_enter_file))
    # optimize the fit i.e. minimize the maximum likelihood estimate.
    bo.estimate_best_kernel_parameters(theta_bounds=[[0.01, 10]])
    logger.info('First kriging -- MLE optimization done')

    # visualize
    visualize_fit(bo)

    # create the new fit
    evaluator_krig = FunctionEvaluator(lambda x: float(bo.evaluate_surrogate_at(x)[0]))
    krig_2 = BayesOpt(cost_function=evaluator_krig, n_dim=n_dim, n_opt=1, n_init=4,
                      kern_function='matern_52', init_strategy=None, u_bound=u_bound,
                      l_bound=l_bound, acq_func='LCB', kappa_strategy=lambda _: 1e9,
                      if_restart=False, do_init=True, acq_func_optimizer=my_acq_optimizer_brute,
                      constraints=[custom_sum_constraint])

    term_criteria = False
    curr_iter = 0
    while not term_criteria:
        # update 10 iteration
        krig_2.update_iter(1)
        if not curr_iter % 5:
            # estimate and set kernel parameters
            krig_2.estimate_best_kernel_parameters(theta_bounds=[[0.01, 10]])
        if not curr_iter % 3:
            # visualize
            visualize_fit(krig_2)
        # calculate estimates of variance -- termination criteria
        term_criteria = (krig_2.curr_points > 25)            # terminate when the current iteration is more than 20
        curr_iter += 1

    # visualize
    visualize_fit(krig_2)
    x_krig_2, y_krig_2 = krig_2.get_total_population()
    io.savemat('/'.join(manual_data_to_enter_file.split('/')[0:-1]) + '/results/' +
               manual_data_to_enter_file.split('/')[-1][0:-4] + '/data.mat', mdict={'x': x_krig_2, 'y': y_krig_2})
    logger.info('Second kriging surface constructed')
