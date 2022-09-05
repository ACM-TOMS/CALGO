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
"""
Example of running a simple bayesian optimization using the class PARyOpt
"""
import logging

import time

from PARyOpt import BayesOpt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from PARyOpt.evaluators import FunctionEvaluator


def cost_function(x: np.array) -> float:
    """
    double exponential function
    :param x: location (single x)
    :return: value
    """
    f = float(np.exp(-(x - 2.) ** 2) + np.exp(-(x - 6.) ** 2 / 10.) + 1.0 / (x ** 2 + 1) + np.exp(-(x - 8.) ** 2) / 12.)
    return -1.0 * f


evaluator = FunctionEvaluator(cost_function)


def visualize_bo(b_opt: BayesOpt) -> None:
    """
    Visualize 1D cost function, surrogate, acquisition function and penalty for inputted bayesian
    optimization class
    :param b_opt: Bayesian optimization class
    :return:
    """
    xs = np.linspace(start=b_opt.l_bound, stop=b_opt.u_bound, num=100, retstep=False)
    _target = np.empty(xs.shape[0])
    _mean = np.empty(xs.shape[0])
    _variance = np.empty(xs.shape[0])
    _acq = np.empty(xs.shape[0])
    _penalty = np.empty(xs.shape[0])
    _s_opt = np.empty(xs.shape[0])
    kappa = b_opt.get_kappa[0](b_opt.curr_iter+1)
    for _i in range(xs.shape[0]):
        _target[_i] = cost_function(xs[_i])
        _mean[_i], var = b_opt.evaluate_surrogate_at(np.asarray([xs[_i]]))
        _variance[_i] = np.sqrt(abs(var))
        _acq[_i] = b_opt.acquisition_function[0](_mean[_i], _variance[_i])
        _penalty[_i] = b_opt.penalty_function(np.asarray([xs[_i]]), kappa)
        _s_opt[_i] = b_opt.acquisition(np.asarray([xs[_i]]), kappa)

    min_s = np.argmin(_s_opt)
    gs = gridspec.GridSpec(4, 1, height_ratios=[4, 1, 1, 1])
    ax0 = plt.subplot(gs[0])
    ax0.plot(xs, _target, 'r-', linewidth=4, label='Actual')
    ax0.plot(xs, _mean, 'g-', linewidth=2, label='Surrogate mean')
    ax0.fill_between(xs, _mean+_variance, _mean-_variance, alpha=0.15)
    ax0.plot(b_opt.next_query_points[0], cost_function(b_opt.next_query_points[0]), 'r*', markersize=14)
    ax0.plot(xs, _mean + _variance, 'b--', label='UCB', alpha=0.2)
    ax0.plot(xs, _mean - _variance, 'k--', label='LCB', alpha=0.2)
    ax0.plot(xs[min_s], _mean[min_s], 'go', markersize=14)
    ax0.legend(loc='upper left')
    ax1 = plt.subplot(gs[1])
    ax1.plot(xs, _acq, 'b-', linewidth=2, label='Acquisition function')
    ax1.plot(xs[min_s], _acq[min_s], 'go', markersize=14)
    ax1.legend(loc='upper left')
    ax2 = plt.subplot(gs[2])
    ax2.plot(xs, _penalty, 'b-', linewidth=2, label='Penalty')
    ax2.plot(xs[min_s], _penalty[min_s], 'go', markersize=14)
    kappa = b_opt.get_kappa[0](b_opt.curr_iter)
    ax2.plot(b_opt.next_query_points[0], b_opt.penalty_function(b_opt.next_query_points[0], kappa=kappa),
             'r*', markersize=14)
    ax2.legend(loc='upper left')
    ax3 = plt.subplot(gs[3])
    ax3.plot(xs, _s_opt, 'b-', linewidth=2, label='Acquisition + Penalty ')
    ax3.plot(b_opt.next_query_points[0], b_opt.acquisition(b_opt.next_query_points[0], kappa), 'r*', markersize=14)
    ax3.plot(xs[min_s], _s_opt[min_s], 'go', markersize=14)
    ax3.legend(loc='upper left')
    mng = plt.get_current_fig_manager()
    # mng.resize(*mng.window.maxsize())
    plt.show()


def user_defined_kappa(curr_iter: int) -> float:
    """
    Sinusoidal variation of kappa
    :param curr_iter: current iteration
    :return: value of kappa
    """
    kappa = 2 * np.exp(-curr_iter * 0.8)
    # kappa = 2*(np.sin(curr_iter * np.pi / 10.) + 1.5)
    # kappa = 2.0
    return kappa


if __name__ == "__main__":

    # logging setup
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler('bayes_opt_log_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S")), mode='a')
    print('bayes_opt_log_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S")))
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info('Performing a simple bayesian optimization')
    # Plotting the target function:
    l_bound = np.asarray([-2.0])
    u_bound = np.asarray([10.0])
    n_dim = 1
    bayesOpt = BayesOpt(cost_function=evaluator, n_dim=n_dim, n_opt=1,
                        n_init=2, kern_function='sqr_exp', init_strategy=None, u_bound=u_bound,
                        l_bound=l_bound, acq_func='LCB', kappa_strategy=user_defined_kappa,
                        if_restart=False)
    logger.info('BO initialized')

    a = np.asarray([-1])
    mean, variance = bayesOpt.evaluate_surrogate_at(a)
    plt.plot(a, mean, '*')
    # sample optimization using scipy:
    # res = minimize(bayesOpt.acquisition, x0=0.5*(l_bound+u_bound), method='Nelder-Mead',
    #                options={'xtol': 1e-4, 'disp': True}, bounds=[(l_bound, u_bound)])

    # a = res.x
    # plt.plot(a, cost_function(a))
    # plt.show()
    for _ in range(20):
        # find next point using surrogate optimization
        tot_pop, tot_func = bayesOpt.get_total_population()
        plt.plot(tot_pop, tot_func, 'yo')

        # visualization of constructed surrogate
        visualize_bo(bayesOpt)
        bayesOpt.update_iter()
        plt.show()

    # plot the locations where evaluation is already performed:
    tot_pop, tot_func = bayesOpt.get_total_population()
    plt.plot(tot_pop, tot_func, 'g*')
    plt.show()
    logger.info('Bayesian optimization completed!')
