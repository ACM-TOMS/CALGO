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
Tests for PARyOpt

This script tests PARyOpt without testing function evaluator, asynchronocity and acquisition functions
"""

from PARyOpt import BayesOpt
import numpy as np
from PARyOpt.evaluators import FunctionEvaluator


def cost_function(x: np.array) -> float:
    """
    double exponential function
    :param x: location (single x)
    :return: value
    """
    f = float(np.exp(-(x - 2.) ** 2) + np.exp(-(x - 6.) ** 2 / 10.) + 1.0 / (x ** 2 + 1) + np.exp(-(x - 8.) ** 2) / 12.)
    return -1.0 * f


def user_defined_kappa(curr_iter: int) -> float:
    """
    Sinusoidal variation of kappa
    :param curr_iter: current iteration
    :return: value of kappa
    """
    kappa = 2.0
    return kappa


def get_next_point(iter_number):
    """
    Manually generate points within the bounds to be fed into BayesOpt instance
    These are uniformly generated random points
    :param iter_number:
    :return:
    """
    next_point = -2.0 + iter_number/50.0 * (10. + 2.)
    return np.asarray([next_point * 0.5])


def init_paryopt(max_points: int, max_failed:int = 0) -> BayesOpt:
    """
    Initialize and return bayes opt with specified max points
    :param max_points:
    :param max_failed:
    :return:
    """
    if max_points <= 2:
        max_points = 2

    evaluator = FunctionEvaluator(cost_function)
    l_bound = np.asarray([-2.0])
    u_bound = np.asarray([10.0])
    n_dim = 1
    bayesOpt = BayesOpt(cost_function=evaluator, n_dim=n_dim, n_opt=1,
                        n_init=2, kern_function='sqr_exp', init_strategy=None, u_bound=u_bound,
                        l_bound=l_bound, acq_func='LCB', kappa_strategy=user_defined_kappa, do_init=False,
                        if_restart=False)
    bayesOpt.set_new_kernel_parameters(theta=0.2)
    bayesOpt.add_point(x=get_next_point(5), y=cost_function(get_next_point(5)), if_check_nearness=True)
    bayesOpt.add_point(x=get_next_point(45), y=cost_function(get_next_point(45)), if_check_nearness=True)
    bayesOpt.update_surrogate()
    iter_num = 0
    # success points
    while bayesOpt.curr_points < max_points:
        next_pt = get_next_point(iter_num + 5)
        bayesOpt.add_point(x=next_pt, y=cost_function(next_pt),
                           if_check_nearness=True)
        bayesOpt.update_surrogate()
        iter_num += 1

    # failed points
    while bayesOpt.curr_points < max_points + max_failed:
        next_pt = get_next_point(iter_num + 5)
        bayesOpt.add_point(x=next_pt, is_failed=True,
                           if_check_nearness=True)
        bayesOpt.update_surrogate()
        iter_num += 1
    return bayesOpt


def test_init() -> None:
    """
    Tests whether the initialization is happening correctly
    :return:
    """
    bo = init_paryopt(2, 0)
    added_population = np.asarray([get_next_point(45), get_next_point(5)])
    # check correct number of points are in the population
    assert bo.curr_points == 2
    assert np.allclose(np.asarray(bo.total_population), added_population)
    # check variance is zero at added points
    assert np.isclose(bo.evaluate_surrogate_at(get_next_point(5))[1], 0.0)
    assert np.isclose(bo.evaluate_surrogate_at(get_next_point(45))[1], 0.0)


def test_success_points() -> None:
    """
    Tests whether the surrogate is interpolating the population
    :return:
    """
    bo = init_paryopt(20, 0)
    population, func = bo.get_total_population()
    for pop in population:
        # check mean and variance
        assert np.isclose(bo.evaluate_surrogate_at(pop)[1], 0.0)
        assert np.isclose(bo.evaluate_surrogate_at(pop)[0], cost_function(pop))


def test_MLE() -> None:
    """
    Tests whether MLE optimization is performing correctly. It should be 0.001 for the set parameters
    Also for failed points
    :return:
    """

    # For all successful surrogate
    bo = init_paryopt(20, 0)
    bo.estimate_best_kernel_parameters(theta_bounds=[[0.0010, 0.01]])
    population, func = bo.get_total_population()
    assert bo.kernel.theta == 0.001
    for pop in population:
        # check mean and variance
        assert np.isclose(bo.evaluate_surrogate_at(pop)[1], 0.0)
        assert np.isclose(bo.evaluate_surrogate_at(pop)[0], cost_function(pop))

    # for surrogate including failed points
    bo = init_paryopt(10, 5)
    bo.estimate_best_kernel_parameters(theta_bounds=[[0.0010, 0.01]])
    population, func = bo.get_total_population()
    assert bo.kernel.theta == 0.001
    for pop in population:
        # check mean and variance
        assert np.isclose(bo.evaluate_surrogate_at(pop)[1], 0.0)
        assert np.isclose(bo.evaluate_surrogate_at(pop)[0], cost_function(pop))

    for pop in bo.failed_query_points:
        # check variance at the failed points. it should be zero
        assert np.isclose(bo.evaluate_surrogate_at(pop, include_failed=True)[1], 0.0)

    for pop in population:
        # the addition of failed points should not affect the interpolation
        assert np.isclose(bo.evaluate_surrogate_at(pop, include_failed=False)[0],
                          bo.evaluate_surrogate_at(pop, include_failed=True)[0])


def test_failed() -> None:
    """
    Tests failed points
    :return:
    """
    bo = init_paryopt(20, 4)
    population, func = bo.get_total_population()
    for pop in population:
        # check mean and variance
        assert np.isclose(bo.evaluate_surrogate_at(pop)[1], 0.0)
        assert np.isclose(bo.evaluate_surrogate_at(pop)[0], cost_function(pop))

    for pop in bo.failed_query_points:
        # check variance at the failed points. it should be zero
        assert np.isclose(bo.evaluate_surrogate_at(pop, include_failed=True)[1], 0.0)

    for pop in population:
        # the addition of failed points should not affect the interpolation
        assert np.isclose(bo.evaluate_surrogate_at(pop, include_failed=False)[0],
                          bo.evaluate_surrogate_at(pop, include_failed=True)[0])

    # randomly generated locations should also have same value with and without the failed points

    locations_for_testing = np.linspace(bo.l_bound, bo.u_bound, num=20)
    for loc in locations_for_testing:
        pop = np.asarray([loc])
        assert np.isclose(bo.evaluate_surrogate_at(pop, include_failed=False)[0],
                          bo.evaluate_surrogate_at(pop, include_failed=True)[0])


