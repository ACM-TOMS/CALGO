"""
Example of running a simple bayesian optimization using the class PARyOpt
"""
import getpass
from typing import List, Any

from sys import platform as _platform

import numpy as np

from PARyOpt import BayesOpt
from PARyOpt.evaluators import FunctionEvaluator, AsyncLocalEvaluator, AsyncSbatchEvaluator
import sys
import os
import io
import logging
import time
if _platform == "darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import gridspec, cm
from mpl_toolkits.mplot3d import axes3d

# for visualization
from PARyOpt.evaluators.connection import Connection


def cost_function(x: np.array) -> float:
    """
    double exponential function
    :param x: location (single x)
    :return: value
    """
    f = 1.0

    for i in range(x.shape[0]):
        coord = x[i]
        f *= np.exp(-(coord - 2.) ** 2) + np.exp(-(coord - 6.) ** 2 / 10.) + 1.0 / (coord ** 2 + 1)

    return -1.0 * f


# for AsyncLocalFunctionEvaluator
def folder_generator(directory, x) -> None:
    # write config file
    with open(os.path.join(directory, 'config.txt'), 'w') as f:
        pass  # write file
    pass


def run_cmd(directory, x) -> List[Any]:
    """
    Command to run on local machine to get the value of cost function at x, in directory
    :param directory:
    :param x:
    :return:
    """
    eval_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'evaluator.py')
    return [sys.executable, eval_path] + list(x[:])


def result_parser(directory, x) -> float:
    """
    Parses the result from a file and returns the cost function. The file is written be the actual cost function
    :param directory:
    :param x:
    :return:
    """
    with open(os.path.join(directory, 'result.txt'), 'r') as f:
        return float(f.readline())


# for AsyncSbatchEvaluator
def sbatch_job_script_gen(lcl_dir: str, x: np.array):
    """
    Job script that needs to be submitted using command sbatch. The whole script should be returned as a string
    :param lcl_dir:
    :param x:
    :return:
    """
    return """#!/bin/bash
#SBATCH --job-name='paryopt_test'
#SBATCH --output='output.txt'
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:02:00

module load python
python evaluator.py {}
""".format(' '.join([str(s) for s in x]))


def result_parser_cluster(local_dir: str, remote_dir: str, conn: Connection, x: np.array) -> float:
    """
    Parse results from a file on a cluster. The file will be located in remote_dir, which is an extension \
    of the files in local_dir.
    :param local_dir:
    :param remote_dir:
    :param conn:
    :param x:
    :return:
    """
    filename = 'myresults.txt'
    fl = io.BytesIO()
    conn.sftp().getfo(os.path.join(remote_dir, filename), fl)
    return float(fl.getvalue().decode(encoding='utf-8').strip())


def visualize_bo(b_opt: BayesOpt) -> None:
    """
    Visualize 1D cost function, surrogate, acquisition function and penalty for inputted bayesian
    optimization class
    :param b_opt: Bayesian optimization class
    :return:
    """
    if b_opt.n_dim == 1 :
        xs = np.linspace(start=b_opt.l_bound, stop=b_opt.u_bound, num=100, retstep=False)
        _target = np.empty(xs.shape[0])
        _mean = np.empty(xs.shape[0])
        _variance = np.empty(xs.shape[0])
        _acq = np.empty(xs.shape[0])
        _penalty = np.empty(xs.shape[0])
        _s_opt = np.empty(xs.shape[0])
        for _i in range(xs.shape[0]):
            _target[_i] = cost_function(xs[_i])
            _mean[_i], var = b_opt.evaluate_surrogate_at(np.asarray([xs[_i]]))
            _variance[_i] = np.sqrt(abs(var))
            _acq[_i] = b_opt.acquisition_function[0](_mean[_i], _variance[_i])
            _penalty[_i] = b_opt.penalty_function(xs[_i], b_opt.get_kappa[0](b_opt.curr_iter+1))
            _s_opt[_i] = b_opt.acquisition(xs[_i])

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
        ax2.plot(b_opt.next_query_points[0], b_opt.penalty_function(b_opt.next_query_points[0], kappa=b_opt.kappa_vec[0]),
                 'r*', markersize=14)
        ax2.legend(loc='upper left')
        ax3 = plt.subplot(gs[3])
        ax3.plot(xs, _s_opt, 'b-', linewidth=2, label='Acquisition + Penalty ')
        ax3.plot(b_opt.next_query_points[0], b_opt.acquisition(b_opt.next_query_points[0]), 'r*', markersize=14)
        ax3.plot(xs[min_s], _s_opt[min_s], 'go', markersize=14)
        ax3.legend(loc='upper left')
        figManager = plt.get_current_fig_manager()
        if _platform != "darwin":
            figManager.window.showMaximized()
        plt.show()
    elif b_opt.n_dim == 2:
        print('2D visualization!')
        xs = np.linspace(start=b_opt.l_bound[0], stop=b_opt.u_bound[0], num=50, retstep=False)
        ys = np.linspace(start=b_opt.l_bound[1], stop=b_opt.u_bound[1], num=50, retstep=False)
        _target = np.empty((xs.shape[0], ys.shape[0]))
        _mean = np.empty((xs.shape[0], ys.shape[0]))
        _variance = np.empty((xs.shape[0], ys.shape[0]))
        _acq = np.empty((xs.shape[0], ys.shape[0]))
        _penalty = np.empty((xs.shape[0], ys.shape[0]))
        _s_opt = np.empty((xs.shape[0], ys.shape[0]))
        X, Y = np.meshgrid(xs, ys)
        kappa = b_opt.get_kappa[0](b_opt.curr_iter)
        for _i, _j in [(x, y) for x in range(xs.shape[0]) for y in range(ys.shape[0])]:
            location = np.asarray([xs[_i], ys[_j]])
            _target[_i, _j] = cost_function(location)
            _mean[_i, _j], var = b_opt.evaluate_surrogate_at(location)
            _variance[_i, _j] = np.sqrt(abs(var))
            _acq[_i, _j] = b_opt.acquisition_function[0](_mean[_i, _j], _variance[_i, _j])
            _penalty[_i, _j] = b_opt.penalty_function(location, kappa)
            _s_opt[_i, _j] = b_opt.acquisition(location, kappa)
        tot_pop, func_pop = b_opt.get_total_population()
        vmin = np.min(np.min(_target))
        vmax = np.max(np.max(_target))
        plt.figure(1)
        ax1 = plt.subplot(3, 2, 1)
        p1 = plt.contourf(X, Y, _target.transpose(), cmap=cm.binary, vmin=vmin, vmax=vmax)

        xs = [p[0] for p in tot_pop]
        ys = [p[1] for p in tot_pop]

        plt.plot(xs, ys, 'r*')
        ax1.set_aspect('equal')
        ax1.set_title('Target function')
        ax2 = plt.subplot(3, 2, 2)
        p2 = plt.contourf(X, Y, _mean.transpose(), cmap=cm.binary, vmin=vmin, vmax=vmax)
        plt.plot(xs, ys, 'r*')
        ax2.set_aspect('equal')
        ax2.set_title('Mean of estimate')
        ax3 = plt.subplot(3, 2, 3)
        p3 = plt.contourf(X, Y, _variance.transpose(), cmap=cm.binary)
        plt.plot(xs, ys, 'r*')
        ax3.set_aspect('equal')
        ax3.set_title('Variance of estimate')
        ax4 = plt.subplot(3, 2, 4)
        p4 = plt.contourf(X, Y, _penalty.transpose(), cmap=cm.binary)
        plt.plot(xs, ys, 'r*')
        ax4.set_aspect('equal')
        ax4.set_title('Penalty function: kappa = {0:.5f}'.format(b_opt.get_kappa[0](b_opt.curr_iter)))
        ax5 = plt.subplot(3, 1, 3)
        p5 = plt.contourf(X, Y, _s_opt.transpose(), cmap=cm.binary)
        plt.plot(xs, ys, 'r*')
        ax5.set_aspect('equal')
        ax5.set_title('Acquisition function for next location')
        figManager = plt.get_current_fig_manager()
        if hasattr(figManager, 'window'):
            figManager.window.showMaximized()
        plt.show(block=True)
        # plot wireframe plots
        fig2 = plt.figure(2)
        ax1 = fig2.add_subplot(3, 2, 1, projection='3d')
        ax1.plot_wireframe(X, Y, _target.transpose())
        ax1.scatter(xs, ys, func_pop, s=14, c='red')
        ax1.set_aspect('equal')
        ax1.set_title('Target function')
        ax2 = fig2.add_subplot(3, 2, 2, projection='3d')
        ax2.plot_wireframe(X, Y, _mean.transpose())
        ax2.scatter(xs, ys, func_pop, s=14, c='red')
        ax2.set_aspect('equal')
        ax2.set_title('Mean of estimate')
        ax3 = fig2.add_subplot(3, 2, 3, projection='3d')
        ax3.plot_wireframe(X, Y, _variance.transpose())
        ax3.set_aspect('equal')
        ax3.set_title('Variance of estimate')
        ax4 = fig2.add_subplot(3, 2, 4, projection='3d')
        ax4.plot_wireframe(X, Y, _penalty.transpose())
        ax4.set_aspect('equal')
        ax4.set_title('Penalty function: kappa = {0:.5f}'.format(b_opt.get_kappa[0](b_opt.curr_iter)))
        ax5 = fig2.add_subplot(3, 1, 3, projection='3d')
        ax5.plot_wireframe(X, Y, _s_opt.transpose())
        ax5.scatter(xs, ys, func_pop, s=14, c='red')
        ax5.set_aspect('equal')
        ax5.set_title('Acquisition function for next location')
        if hasattr(figManager, 'window'):
            figManager.window.showMaximized()
        plt.show(block=True)

    else:
        pass


def user_defined_kappa(curr_iter: int, freq: float = 10., t_const: float = 0.5) -> float:
    """
    kappa is the parameter in acquisition function. Current design of acquisition function
    requires large kappa value for exploration and small kappa value for exploitation
    :param curr_iter: current iteration
    :param freq : frequency of sine input
    :param t_const: time constant for exponential
    :return: value of kappa
    """
    # kappa = 2 * np.exp(curr_iter * 0.8)
    kappa = 40.5 * (np.sin(curr_iter * np.pi / freq) + 1.5) * np.exp(-t_const*curr_iter)
    # kappa = 2.0
    return kappa


def perform_optimization(n_dim: int, async_fraction: float, n_opt: int = 4) -> None:
    """
    Performs bayesian optimization for given dimensionality and async fraction
    :param n_dim: dimension of problem
    :param async_fraction: asynchronous fraction
    :param n_opt: optima per iteration
    :return:
    """
    # logging setup
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler('async_opt_{}_ndim_{}_frac_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S"), n_dim,
                                                                       async_fraction), mode='a')

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info('Performing a simple bayesian optimization')

    l_bound = -4.*np.ones((n_dim,))
    u_bound = 10.*np.ones((n_dim,))
    # NOTE THIS
    n_init = 2 * n_dim
    iter_max = 8

    # parallel, asynchronous, out-of-script
    evaluator = AsyncLocalEvaluator(job_generator=folder_generator,
                                    run_cmd_generator=run_cmd,
                                    parse_result=result_parser,
                                    required_fraction=async_fraction)

    # serial, in-script
    # evaluator = FunctionEvaluator(cost_function)

    # parallel, asynchronous, on Comet
    # evaluator = AsyncSbatchEvaluator(Host(getpass.getuser(), 'comet.sdsc.edu'), folder_generator,
    #                                  sbatch_job_script_gen, lcl_parse_result=VALUE_FROM_FILE('result.txt'))

    my_kappa_funcs = []
    for j in range(n_opt):
        my_kappa_funcs.append(lambda curr_iter, freq=10.*(j*j+2), t_const=0.8/(1. + j):
                              user_defined_kappa(curr_iter, freq=freq, t_const=t_const))

    bayesOpt = BayesOpt(cost_function=evaluator, n_dim=n_dim, n_opt=n_opt,
                        n_init=n_init, kern_function='matern_52', init_strategy=None, u_bound=u_bound,
                        l_bound=l_bound, acq_func='LCB', kappa_strategy=my_kappa_funcs,
                        if_restart=False)
    logger.info('BO initialized')

    # a = np.asarray([-1])
    # mean, variance = bayesOpt.evaluate_surrogate_at(a)
    # plt.plot(a, mean, '*')
    # sample optimization using scipy:
    # res = minimize(bayesOpt.acquisition, x0=0.5*(l_bound+u_bound), method='Nelder-Mead',
    #                options={'xtol': 1e-4, 'disp': True}, bounds=[(l_bound, u_bound)])

    # a = res.x
    # plt.plot(a, cost_function(a))
    # plt.show()
    for _ in range(iter_max):
        bayesOpt.update_iter()
        bayesOpt.estimate_best_kernel_parameters(theta_bounds=[[0.05, 10.]])
        visualize_bo(bayesOpt)

    # export cost function evaluations to a CSV file
    bayesOpt.export_csv('data.csv')
    bayesOpt.estimate_best_kernel_parameters(theta_bounds=[[0.05, 10.]])

    # visualization of constructed surrogate
    visualize_bo(bayesOpt)

    logger.info('Bayesian optimization completed!')
    logger.removeHandler(fh)


if __name__ == "__main__":
    perform_optimization(n_dim=2, async_fraction=0.5, n_opt=4)
    perform_optimization(n_dim=2, async_fraction=1.0, n_opt=8)
    perform_optimization(n_dim=4, async_fraction=0.5, n_opt=8)
    perform_optimization(n_dim=4, async_fraction=1.0, n_opt=16)
    perform_optimization(n_dim=10, async_fraction=1.0, n_opt=8)
