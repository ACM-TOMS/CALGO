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
Main class for bayesian optimization
"""
import logging
import math
import os
import pickle
from scipy import io
import traceback
import warnings
import csv  # for data export
from typing import Callable, Union, List, Tuple

import numpy as np
from scipy import optimize

import PARyOpt.acquisition_functions as acqf
import PARyOpt.kernel as kf
import PARyOpt.utils as ut


# logging
logger = logging.getLogger(__name__)


def _get_kernel_function(name: str) -> kf.KernelFunction:
    """
    Returns a kernel function corresponding to a name.
    Raises a ValueError if name is unknown.

    :param name: name of kernel function
    :return: function corresponding to name
    """
    if name == 'sqr_exp':
        return kf.SquaredExponential()
    elif name == 'matern_52':
        return kf.Matern52()
    elif name == 'matern_32':
        return kf.Matern32()
    else:
        raise ValueError('Unknown / unimplemented kernel function')


def _get_acquisition_function(name: str) -> Callable[[float, float, float, float], float]:
    """
    Returns an acquisition function corresponding to a name.
    Raises a ValueError if name is unknown.

    :param name:
    :return:
    """
    if name.upper() == 'LCB':
        # Lower confidence bound : for minimization problems
        return acqf.lower_confidence_bound
    elif name.upper() == 'UCB':
        # Upper confidence bound : for maximization problems
        return acqf.upper_confidence_bound
    elif name.upper() == 'PI':
        # Probability of improvement
        return acqf.probability_improvement
    elif name.upper() == 'EI':
        # expected improvement
        return acqf.expected_improvement
    else:
        raise ValueError('Unknown / unimplemented acquisition function')


def _resolve_acquisition_list(acq_list: List[Union[str, Callable]]= None) -> List[Callable]:
    """
    Resolves the list of acquisition functions and returns a list of acquisition function

    :param acq_list: list of strings of functions
    :return: list of functions
    """
    acq_function_list = []
    for acq_func in acq_list:
        if type(acq_func) is str:
            acq_function_list += [_get_acquisition_function(acq_func)]
        elif callable(acq_func):
            acq_function_list += [acq_func]
        else:
            raise ValueError('Unknown / unimplemented acquisition function')
    return acq_function_list


def _constant_kappa(curr_iter: int) -> float:
    """
    Return a constant kappa
    This should also be a sample of kappa strategy

    :param curr_iter: current iteration
    :return:
    """
    assert curr_iter >= 0
    return 2.0


def _init_strategy_lhs(n_dim: int, n_init: int, l_bound: np.array, u_bound: np.array) -> List[np.array]:
    """
    Generate initial locations for the population.
    The locations are sampled between 0 and 1 using Latin Hypercube Sampling, then scaled to the bounds.

    :param n_init: number of initial points
    :param n_dim: dimensionality
    :return: an array of shape (n_init, n_dim)
    """
    normalized_population = ut.lhs(n_dim, samples=n_init, criterion='cm', iterations=20)  # type: np.ndarray
    real_population = l_bound + normalized_population * (u_bound - l_bound)
    real_population = [p for p in real_population]
    return real_population


def _default_acq_optimizer(func: Callable[[np.array], float], x0: np.array, l_bound: np.array, u_bound: np.array) -> \
        np.array:
    """
    Default acquisition function optimizer : Uses the Powell method in scipy.optimize

    :param func: acquisition function
    :param x0: initial guess
    :param l_bound : lower bound
    :param u_bound : upper bound
    :return: Optimum location (array)
    """
    # Another example using the brute search
    # res = optimize.brute(func, bounds, Ns=50, full_output=True)
    # return res[0]

    bounds = list(zip(l_bound, u_bound))
    res = optimize.minimize(fun=func, x0=x0, bounds=bounds, method='Powell')
    return res.x


def _invert_matrix_mult(A: np.array, b: np.array) -> Tuple[np.array, np.array]:
    """
    Invert the input matrix and multiply with b to get coeff, i.e., returns inv(A) and A_inv_b
    :param A: Matrix
    :param b: vector
    :return: A^(-1) and A^(-1) b
    """
    n = A.shape[0]
    assert A.shape[0] == b.shape[0]
    try:
        # cholesky decomposition
        logger.debug('Eigen values of K are: {}'.format(np.linalg.eigvals(A)))
        l_cholesky = np.linalg.cholesky(A)
        l_inv = np.linalg.solve(l_cholesky, np.identity(n))
        A_inv = l_inv.transpose().dot(l_inv)
    except np.linalg.linalg.LinAlgError:
        logger.debug('Using direct inverse and not Cholesky decomposition. near singular matrix')
        A_inv = np.linalg.inv(A)

    A_inv_b = A_inv.dot(b)

    return A_inv, A_inv_b


class BayesOpt:
    """
    Bayesian optimization class.

    :var curr_iter: iteration number of optimizer
    :var n_surrogate_opt: number of optima to get during surrogate optimization
    :var n_init: initial population size (>2)
    :var n_dim: dimensions of optimization
    :var l_bound, u_bound: lower and upper bounds on variables
    :var total_population: total set of population: list of arrays
    :var func_vector: functional values for the population : scalar for each population
    :var K: covariance matrix
    :var K_inv: inverse of covariance matrix
    :var K_inv_y: K_inv * func_vector : pre-computation to save costly / unstable inversion
    :var acquisition_function: acquisition function to optimize the surrogate
    :var cost_function: cost function to MINIMIZE: should be able to take a vector of locations
    :var constraints: list of constraint functions (only points that satisfy these functions will be visited) \
    NOT enforced in the default init strategy!!
    :var if_restart: whether it should restart or not
    :var restart_filename: file to restart from. it will be `opt_state.dat` if nothing is specified
    :var acq_func_optimizer: optimizer for acquisition function. Should take in function, initial guess, bounds and \
    function derivative. Returns the optimal location
    """

    def __init__(self,
                 cost_function: Union[
                     Callable[[List[np.array]],
                              Tuple[List[Tuple[np.array, float]], List[np.array], List[np.array]]],
                     Callable[[List[np.array], List[np.array]],
                              Tuple[List[Tuple[np.array, float]], List[np.array], List[np.array]]]],
                 l_bound: np.array, u_bound: np.array, n_dim: int,
                 n_opt: int = 1, n_init: int = 0,
                 init_strategy: Callable[[int, int, np.array, np.array], List[np.array]] = None,
                 do_init: bool = True,
                 kern_function: Union[str, kf.KernelFunction] = None,
                 acq_func: Union[str, Callable[[float, float, float, float], float]] = None,
                 acq_func_optimizer: Callable = None,
                 kappa_strategy: Union[List[Callable[[int], float]], Callable[[int], float], None] = None,
                 constraints: List[Callable[[np.array], bool]] = list(),
                 if_restart: bool = None, restart_filename: str = 'opt_state.dat') -> None:

        assert l_bound.shape == u_bound.shape, "The supplied bounds have unequal shapes"
        assert len(l_bound) == n_dim, "The bounds have different dimensionality than n_dim"
        assert n_opt >= 1, "Asking for less than 1 optima per iteration, NOT possible. Terminating!"
        if do_init:
            assert n_init > 1, "Initial points are less than 2, not possible to build a surrogate"

        self.n_surrogate_opt = n_opt
        self.n_dim = n_dim
        self.n_init = n_init
        self.u_bound = u_bound
        self.l_bound = l_bound

        assert kern_function is not None, "No kernel function specified"
        if type(kern_function) is str:
            self.kernel = _get_kernel_function(kern_function)
        else:
            assert hasattr(kern_function, "eval"), "Kernel function does not implement the KernelFunction interface"
            self.kernel = kern_function

        # resolve acquisition functions
        assert callable(acq_func) or type(acq_func) is str or type(acq_func) is list, "Supplied acquisition function " \
                                                                                      "is neither callable nor a " \
                                                                                      "string" \
                                                                                      "nor a list of callables/string"
        if type(acq_func) is str or callable(acq_func):
            # only one string/function is provided, make it into a list of strings/functions
            acq_func_list = [acq_func for _ in range(self.n_surrogate_opt)]
            self.acquisition_function = _resolve_acquisition_list(acq_func_list)
        elif type(acq_func) is list:
            assert len(acq_func) == self.n_surrogate_opt, "Number of acquisition functions supplied do not " \
                                                          "match the requested optima per iteration"
            self.acquisition_function = _resolve_acquisition_list(acq_func)
        else:
            raise ValueError("Problem with acquisition function")

        # acquisition function optimizer
        assert acq_func_optimizer is None or callable(acq_func_optimizer), \
            "Supplied acquition function optimizer is not callable"
        self.acq_func_optimizer = acq_func_optimizer if acq_func_optimizer is not None else _default_acq_optimizer

        assert cost_function is not None, "No cost function supplied, terminating!"
        self.cost_function = cost_function

        # managed variables
        # current iteration number
        self.curr_iter = 0  # type: int
        # number of points in the surrogate
        self.curr_points = 0    # type: int
        self.total_population = list()  # type: List[np.array]

        # allocate function values for total population
        self.func_vector = list()  # type: List[float]

        # covariance matrix & inverse
        self.K = np.identity(self.curr_points)
        self.K_inv = np.identity(self.curr_points)
        self.K_success = np.identity(self.curr_points)
        self.K_inv_success = np.identity(self.curr_points)

        # inverse times the function vector :
        self.K_inv_y = np.empty((self.curr_points, 1))
        self.K_inv_y_success = np.empty((self.curr_points, 1))

        # restart
        self.restart_file_path = restart_filename
        if if_restart is None:
            if_restart = os.path.exists(self.restart_file_path)
        self.if_restart = if_restart
        self.inside_iter = False

        # perform initialization based on strategy provided:
        if do_init and init_strategy is None:
            self.init_strategy = _init_strategy_lhs
            if len(constraints) > 0:
                warnings.warn("User-specified constraints are not respected by the default init_strategy. "
                              "Initialization will probably fail.", category=RuntimeWarning)
        else:
            self.init_strategy = init_strategy

        # kappa strategy:
        # if not given, default to _constant_kappa
        if kappa_strategy is None:
            kappa_strategy = _constant_kappa

        if type(kappa_strategy) is list:
            # if the user passed a list, they have specified a kappa strategy for each surrogate optima
            # verify the number of functions matches and use it directly
            if not len(kappa_strategy) == self.n_surrogate_opt:
                raise ValueError("Number of kappa strategies does not match the number of surrogate optima")
            self.get_kappa = kappa_strategy
        elif callable(kappa_strategy):
            # kappa_strategy is a function: repeat it n_surrogate_opt times
            self.get_kappa = [kappa_strategy for _ in range(self.n_surrogate_opt)]
        else:
            raise ValueError("Invalid kappa_strategy (expected callable or list of callables)")

        # constraints
        self.constraints = constraints

        # initialize next query points and pending query points
        self.next_query_points = [np.zeros(self.n_dim) for _ in range(self.n_surrogate_opt)]
        self.pending_query_points = list()
        self.failed_query_points = list()
        self.failed_evaluations = list()

        ########### logging ############
        logger.info('dimensionality: {}'.format(n_dim))
        logger.info('Bounds: l_bound: {}, u_bound: {}'.format(l_bound, u_bound))
        logger.info('kernel function: {}'.format(kern_function if type(kern_function) is str
                                                 else type(kern_function).__name__))
        logger.info('acquisition function: {}'. format(acq_func if type(acq_func) is str else acq_func.__name__))
        logger.info('kernel function : theta: {}, theta0: {}'.format(self.kernel.theta, self.kernel.theta0))
        logger.info('Init finished')

        if do_init:
            self._initialize()

        assert all(u_bound >= l_bound), "all the upper bound values are not greater than lower bound values"

    def _initialize(self) -> None:
        """
        Initialize the population based on a random number, between the upper and lower bound

        :return: None
        """
        if self.if_restart:
            self._restart_from_file()
        else:
            # randomly initialize population
            population = self.init_strategy(self.n_dim, self.n_init, self.l_bound, self.u_bound)

            # save pending and failed directly (since this is initialization)
            completed = self._require_all_cost_function(population)
            for x, y in completed:
                assert len(x.shape) == 1 and x.shape[0] == self.n_dim
                self.add_point(x, y, if_check_nearness=True)

            self.update_surrogate()
            # note that np.dot(self.K, self.K_inv) = np.identity
        if self.curr_points < 2:
            raise RuntimeError('Initialization strategy failed to provide at least two points. Cannot continue.')

    def add_point(self, x: np.array, y: float = None, if_check_nearness: bool = True, is_failed: bool = False) -> None:
        """
        Adds ONLY ONE point to the current set of data we have.
        if_check_nearness specifies if it is needed to check nearness constraints before adding
        the point into total_population. This should be true if the user is
        manually adding points into the population.
        Updates following variables: K, K_inv, total_population, func_vector
        If the point is a failed evaluation, then it does not add into total_population and func_vector

        :param x: array of shape (1, self.n_dim)
        :param y: cost function value at that location
        :param if_check_nearness: boolean, whether to check nearness or not.
        :param is_failed: boolean, whether the point to add is a failed evaluation or not
        :return: None
        """
        assert len(x.shape) == 1, "Adding more than one point using add_point. NOT possible, Terminating!"
        assert x.shape[0] == self.n_dim
        if (not if_check_nearness) or (if_check_nearness and not self._is_near_to_population(x)):
            self.curr_points += 1
            self._add_point_to_population(x, is_failed)
            k_vec = self._calc_k_vec_at(x, include_failed=True)
            if not is_failed:
                if y is None:
                    logger.info('add_point(..) evaluating new point, no cost function value supplied')
                    completed = self._require_all_cost_function([x])
                    y = completed[0]
                # add new location to total population
                self._add_point_to_func_vector(y)
            self._add_row_to_covariance(k_vec, is_failed)
            logger.debug("New point added {}!. Current size of population: {}".format(
                'to population' if not is_failed else 'from failed list', self.curr_points))
        else:
            logger.info('The added point ({}) is close to other existing points. Not adding into total population!'.
                        format(x))

    def _calc_log_likelihood(self) -> float:
        """
        Calculates the log-likelihood for the current kernel shape parameters (Maximum Likelihood Estimate, or MLE).
        This operation is done ONLY on self.K_success, i.e., for points with successful function evaluation.

        :return: log-likelihood value for current kernel shape parameters
        """
        n = 1.0 * self.K_success.shape[0]
        '''
        sum_log_eig = 0.0
        eig_val, _ = np.linalg.eig(self.K)         
        for i in range(0, n):
            sum_log_eig += log(eig_val[i])

        sum_log_eig /= n
        '''

        det_K = np.linalg.det(self.K_success)

        if det_K <= 1e-14:
            logger.warning('Negative or zero det(K) for theta = {}, returning penalty'.format(self.kernel.theta))
            return 1.e16

        # from 'Mongillo, Michael. "Choosing basis functions and shape parameters for radial basis function methods."
        # SIAM Undergraduate Research Online 4 (2011): 190-209.'
        sum_log_eig = math.log(det_K) / n
        yT_times_K_inv = np.dot(self.func_vector, self.K_inv_y_success)
        log_yT_c = math.log(yT_times_K_inv)
        mle = log_yT_c + sum_log_eig
        logger.debug("MLE = {}, for theta = {}".format(mle, self.kernel.theta))
        return mle

    def estimate_best_kernel_parameters(self, theta_bounds=None) -> None:
        """
        Calculates and *sets* the best shape-parameter/characteristic length-scale for RBF kernel function and
        applies it to the current model.

        :param theta_bounds: array of bounds for theta
        :return: None
        """
        pre_existing_points = self.curr_points
        logger.debug('Estimating kernel parameters!, {} points already existing'.format(self.curr_points))
        old_parameters = (self.kernel.theta0, self.kernel.theta)
        if theta_bounds is None:
            theta_bounds = [[0.1, 12.0]]

        def calc_mle(x1: float, x2: float = None) -> float:
            """
            Calculate the Maximum likelihood estimate

            :param x1: value of theta for kernel function
            :param x2: value of theta0, generally 1.0
            :return:
            """
            if x2 is None:
                x2 = 1.0
            self.set_new_kernel_parameters(theta=float(x1), theta0=x2)
            mle = self._calc_log_likelihood()
            return mle
        try:
            theta_opt = optimize.brute(calc_mle, theta_bounds, Ns=30, full_output=False, finish=None)
            self.set_new_kernel_parameters(theta=theta_opt, theta0=1.0)
            logger.info('Kernel function parameter optimized: theta: {}, theta0: {}'.format(self.kernel.theta,
                                                                                            self.kernel.theta0))
        except Exception as _:
            logger.error(traceback.format_exc())
            logger.warning('bayes opt hyper parameter optimization failed! Resetting to old parameters')
            self.set_new_kernel_parameters(theta=old_parameters[1], theta0=old_parameters[0])

        assert pre_existing_points == self.curr_points, 'Estimation of MLE is losing points'

    def update_iter(self, iter_max: int = 1) -> None:
        """
        Finds the next query point and updates the data for the next iteration

        :return: None
        """
        for _ in range(iter_max):
            if not self.inside_iter:
                self.next_query_points = self._find_next_query_point()
                self.inside_iter = True
                self._save_to_file()

            # unpack directly into pending
            completed, self.pending_query_points, failed =\
                self.cost_function(self.next_query_points, self.pending_query_points)

            logger.info('Cost function evaluated! Number completed: {}, failed: {} and pending: {}'.format(
                len(completed), len(self.pending_query_points), len(failed)))

            for x, y in completed:
                self.add_point(x, y, if_check_nearness=True, is_failed=False)

            # add the failed points as interpolated points to covariance matrix
            for x in failed:
                self.add_point(x, if_check_nearness=True, is_failed=True)

            # update covariance matrix inverse after adding all points
            self.update_surrogate()
            self.curr_iter += 1

            self.inside_iter = False
            self._save_to_file()

            if len(completed) == 0:
                msg = 'Failed to evaluate any new points during iteration ' + str(self.curr_iter) + '.'
                logger.warning(msg)

        self._save_to_file()
        logger.info('{} iterations completed, total population: {}, number of iterations: {}!'.
                    format(iter_max, self.curr_points, self.curr_iter))

    def _find_next_query_point(self) -> List[np.array]:
        """

        :return: list x values (length up to self.n_surrogate_opt - may be less)
        """
        next_points = []
        avoid_points = []
        max_tries = 10
        for i in range(self.n_surrogate_opt):
            tries = 0
            while tries < max_tries:
                tries += 1
                kappa = self.get_kappa[i](self.curr_iter)
                logger.info('opt = {0:d}, kappa = {1:.2f}'.format(i, kappa))

                def opt_func(x1: np.array) -> float:
                    """
                    function to optimize (acquisition with other arguments filled in)
                    :param x1:
                    :return:
                    """
                    test_pt = np.asarray([x1]) if len(x1.shape) == 0 else x1
                    return self.acquisition(test_pt, kappa, avoid=next_points+avoid_points, opt_indx=i)

                x = self.acq_func_optimizer(func=opt_func,
                                            x0=np.mean([self.l_bound, self.u_bound], axis=0),
                                            l_bound=self.l_bound, u_bound=self.u_bound)
                x = np.reshape(x, (self.n_dim,))          # for cases of 1 dimension
                # never look for this location, even if it is not accepted within the tries. If it is accepted, it will
                # be included in next_points as well, which will increase the penalty at these locations
                avoid_points.append(x)
                if not self._is_near_to_population(x, population=self.total_population + next_points +
                                                                 self.pending_query_points + self.failed_query_points)\
                        and (all(x <= self.u_bound) and all(x >= self.l_bound))\
                        and all([c(x) for c in self.constraints]):
                    next_points.append(x)
                    logger.info('Next location: ' + str(x) + ', (acq value: ' +
                                str(self.acquisition(x, kappa, avoid=next_points, opt_indx=i)) + ')')
                    break
            if tries >= max_tries:
                logger.warning('Could not find surrogate optima ' + str(i) + ' (all points near population)')
        # find unique set in the list of next_points, this will make its size less than n_surrogate_opt
        a = [tuple(x.tolist()) for x in next_points]
        next_points = [np.array(x) for x in set(a)]
        return next_points

    def acquisition(self, x: np.array, kappa: float, avoid: List[np.array] = None, opt_indx: int = 0) -> float:
        """
        Calculates the acquisition function + penalty function at the given query location.
        If the query point is outside the bounds, it will return an infinity
        This is the function that needs to be optimized.

        :param x: location to evaluate acquisition function : shape: (1, self.n_dim)
        :param kappa: kappa value to pass to acquisition function
        :param avoid: list of x values to avoid (forwarded to penalty function)
        :param opt_indx: index of acquisition function that is being called
        :return: scalar acquisition function at x
        """

        assert len(x.shape) == 1, "acquisition received more than one({}) query point!".format(len(x.shape))
        if avoid is None:
            avoid = list()
        if (x > self.u_bound).any() or (x < self.l_bound).any():
            return 1.e16
        if not all([c(x) for c in self.constraints]):
            return 1.e16

        mean, variance = self.evaluate_surrogate_at(x, include_failed=True)
        curr_best = np.min(self.func_vector)
        return self.acquisition_function[opt_indx](mean, np.sqrt(abs(variance)), curr_best, kappa)\
            + 10.0 * self.penalty_function(x, kappa, avoid)

    def penalty_function(self, x: np.array, kappa: float, avoid: List[np.array] = list()) -> float:
        """
        Calculates the penalty function, to prevent local optima to repeat

        :param x: location
        :param kappa: kappa value used in acquisition function
        :param avoid: extra points to avoid (other than self.total_pop/self.pending/self.failed)
        :return: penalty function value at x
        """
        penalty_val = 0.0
        penalty_wid = max(kappa, 1e-2)
        # penalize existing points and avoid list
        for pop in (self.total_population + self.pending_query_points + self.failed_query_points):
            penalty_val += ut.pdf_normal(x, pop, np.sqrt(penalty_wid))
        for pop in avoid:
            penalty_val += 5.0 * ut.pdf_normal(x, pop, np.sqrt(penalty_wid))
        # penalize points close to boundaries
        # iterate through all the boundaries
        for bndr_dim in range(self.n_dim):
            # penalty at the upper bound
            penalty_val += 5.0 * ut.pdf_normal(x[bndr_dim], self.u_bound[bndr_dim], 1. * np.sqrt(penalty_wid))
            # penalty at the lower bound
            penalty_val += 5.0 * ut.pdf_normal(x[bndr_dim], self.l_bound[bndr_dim], 1. * np.sqrt(penalty_wid))
        return penalty_val

    def _save_to_file(self) -> None:
        """
        saves current data to file

        :return: None
        """
        data = {
            "n_dim": self.n_dim,
            "curr_iter": self.curr_iter,
            "curr_points": self.curr_points,
            "total_population": self.total_population,
            "func_vector": self.func_vector,
            "inside_iter": self.inside_iter,
            "next_query_points": self.next_query_points,
            "pending_query_points": self.pending_query_points,
            "failed_query_points": self.failed_query_points,
            "u_bound": self.u_bound,
            "l_bound": self.l_bound,
            "theta0": self.kernel.theta0,
            "theta": self.kernel.theta,
            "K": self.K,
            "K_inv": self.K_inv,
            "K_inv_y": self.K_inv_y,
        }
        with open(self.restart_file_path, 'wb') as f:
            pickle.dump(data, f)

        io.savemat(self.restart_file_path[:-4] + '.mat', data)

        if self.inside_iter:
            logger.info("Saved optimizer state from iter " + str(self.curr_iter) + " to " + str(self.curr_iter + 1))
        else:
            logger.info("Saved optimizer state for iteration " + str(self.curr_iter) + ".")

    def _restart_from_file(self) -> None:
        """
        Restarts from a file

        :return: None, updates all the variables
        """
        with open(self.restart_file_path, 'rb') as f:
            data = pickle.load(f)

        assert self.n_dim == data["n_dim"], "Dimensionality of restart and initialization are not equal"
        self.curr_iter = data["curr_iter"]                          # type: int
        curr_points = data.get("curr_points")                       # type: float
        total_population = data["total_population"]                 # type: List[np.array]
        func_vector = data["func_vector"]                           # type: List[float]
        self.inside_iter = data["inside_iter"]                      # type: bool
        self.next_query_points = data["next_query_points"]          # type: List[np.array]
        self.pending_query_points = data["pending_query_points"]    # type: List[np.array]
        failed_query_points = data["failed_query_points"]           # type: List[np.array]
        self.u_bound = data["u_bound"]                              # type: List[float]
        self.l_bound = data["l_bound"]                              # type: List[float]
        self.kernel.theta = data["theta"]                           # type: float
        self.kernel.theta0 = data["theta0"]                         # type: float

        # since the population is prepended for successful evaluations, restart should add points
        # in the reverse order
        for pop, func in zip(total_population[::-1], func_vector[::-1]):
            # add point to the population, no need to check nearness
            self.add_point(pop, func, if_check_nearness=False, is_failed=False)
        self._evaluate_surrogate_for_failed_list()
        for pop in failed_query_points:
            self.add_point(pop, if_check_nearness=True, is_failed=True)

        assert self.curr_points == curr_points, "Added points({}) and num of points({}) from restart file " \
                                                "do not match".format(self.curr_points, curr_points)

        self.update_surrogate()
        logger.info("Resumed on iteration " + str(self.curr_iter) + ".")

    def evaluate_surrogate_at(self, x: np.array, include_failed: bool = False) -> tuple:
        """
        Evaluates the surrogate at given point
        Taken from : https://arxiv.org/pdf/1012.2599.pdf :
        "A Tutorial on Bayesian Optimization of Expensive Cost Functions, with Application to
        Active User Modeling and Hierarchical Reinforcement Learning"
        Note that it returns sigma^2 and not sigma

        The mean of the surrogate will not be affected by the failed points, they will only affect
        the variance of the surrogate

        :param x: location to evaluate the surrogate rbf approximation
        :param include_failed: whether to include failed points for surrogate
        :return: mean(mu), variance (sigma^2)
        """
        k_vec = self._calc_k_vec_at(x, include_failed=include_failed)  # row vector
        k_vec_transpose = np.asarray([k_vec])  # column vector

        k_cov_diag = self.kernel.eval(x, x)

        if include_failed:
            K_inv_y = self.K_inv_y
            K_inv = self.K_inv
        else:
            K_inv = self.K_inv_success
            K_inv_y = self.K_inv_y_success

        mean = k_vec_transpose.dot(K_inv_y)
        variance = k_cov_diag - k_vec_transpose.dot(K_inv.dot(k_vec))
        return mean, variance

    def _evaluate_surrogate_for_failed_list(self) -> None:
        """
        Evaluates the surrogate values for the failed list. This is used to include failed
        evaluations into the surrogate.
        The surrogate is updated by calling self.update_surrogate before evaluating the interpolated value.

        Updates the class variable self.failed_evaluations
        :return: None
        """
        self.failed_evaluations = []
        for failed_point in self.failed_query_points:
            mean, _ = self.evaluate_surrogate_at(failed_point, include_failed=False)
            self.failed_evaluations.append(mean)
        logger.debug('Surrogate evaluated for failed points.')

    def _calc_k_vec_at(self, x: np.array, include_failed: bool) -> np.array:
        """
        Calculates the k_vector for any arbitrary point. This function can be useful
        when updating the covariance matrix (per iteration) or when evaluating the
        surrogate function at x

        :param x: location to calculate k_vector
        :return: vector of dimension curr_points
        """
        # if not include_failed:
        #     k_vec = np.asarray([self.kernel.eval(x, y) for y in self.total_population])
        # else:
        #     k_vec = np.asarray([self.kernel.eval(x, y) for y in self.total_population] +
        #                        [self.kernel.eval(x, y) for y in self.failed_query_points])

        # lengths of successful and failed points
        pop_num = len(self.total_population)
        fail_num = len(self.failed_query_points)
        if include_failed:
            k_vec = np.empty(self.curr_points)
        else:
            k_vec = np.empty(pop_num)
        # fill k_vec for successful population
        for pop in range(pop_num):
            k_vec[pop] = self.kernel.eval(x, self.total_population[pop])

        if include_failed:
            # for calculating kernel with failed points, the length scale has to be reduced
            self.kernel.theta /= 10
            for pop in range(fail_num):
                k_vec[pop_num + pop] = self.kernel.eval(x, self.failed_query_points[pop])
            # revert the length scale of the kernel
            self.kernel.theta *= 10
        return k_vec

    def _add_point_to_population(self, x: np.array, is_failed: bool = False) -> None:
        """
        adds an element to the total population
        In the case of a failed location, it appends to self.failed_query_points
        In the case of successful evaluation, it prepends to the self.total_population

        The above methodology is used to make calculations of k_vec easier, especially if failed points have smaller
        correlation length. The order of rows in covariance matrix K can be consistent with the
        "list: self.total_population + self.failed_query_points"

        :rtype: None
        :param is_failed: whether we are adding a failed location or not
        :param x: new location to be added to the total population
        :return: None
        """
        if is_failed:
            self.failed_query_points.append(x)
        else:
            self.total_population.insert(0, x)
#        assert len(self.total_population)+len(self.failed_query_points) == self.curr_points, "total population is " \
#                                                                                          "not the same as current" \
#                                                                                          " number of points"

    def _add_point_to_func_vector(self, y: float) -> None:
        """
        prepends an element to the func_vector. This is to consider the order of elements in self.total_population
        (see self._add_point_to_population() for exact details)

        :param y: value of cost function to be added to func_vector
        :return: None
        """
        self.func_vector.insert(0, float(y))

    def _add_row_to_covariance(self, k_vec_new: np.array, is_failed: bool = False) -> None:
        """
        Adds a row to covariance matrix
        ENSURE that self.curr_points is updated before calling this function

        :param k_vec_new: new row to be added to covariance matrix
        :return: None
        """
        # preallocate updated covariance matrix
        k_new = np.identity(self.curr_points)
        if is_failed:
            # failed point, so add to the bottom of the covariance matrix
            # gather from previous covariance matrix and put in the upper part of new covariance matrix
            k_new[0:self.curr_points - 1, 0:self.curr_points - 1] = self.K.copy()
            # fill in the new row/column
            k_new[-1, :] = k_vec_new
            k_new[:, -1] = k_vec_new
        else:
            # normal point, add to the top of the covariance matrix
            # gather from previous covariance matrix and put in the lower part of new covariance matrix
            k_new[1:self.curr_points, 1:self.curr_points] = self.K.copy()
            # fill in the new row/column
            k_new[0, :] = k_vec_new
            k_new[:, 0] = k_vec_new
        # update the class variable
        self.K = k_new.copy()

    def update_surrogate(self) -> None:
        """
        updates the inverse of covariance function.
        In order to incorporate failed locations as interpolated values, this update is 2 staged:

        1. Create the surrogate (K,K_inv) with only the successfully evaluated points

        2. a) Evaluate the mean of the surrogate at the failed points
           b) Incorporate the interpolated failed locations into the surrogate and update the surrogate (K,K_inv)

        Updates the values of self.K_success, self.K_inv_success, self.K_inv_y_success, self.K_inv and self.K_inv_y

        :return: None
        """
        pop_num = len(self.total_population)
        fail_num = len(self.failed_query_points)
        assert self.K.shape[0] == self.curr_points
        assert self.curr_points == pop_num + fail_num

        # create the surrogate without failed points
        self.K_success = self.K[0:pop_num, 0:pop_num]
        self.K_inv_success, self.K_inv_y_success = _invert_matrix_mult(self.K_success, np.asarray(self.func_vector))
        logger.info('Surrogate without failed points created! ({} points) '.format(pop_num))

        # evaluate the interpolated values for failed points
        self._evaluate_surrogate_for_failed_list()
        # create the surrogate with failed points
        self.K_inv, self.K_inv_y = _invert_matrix_mult(self.K, np.asarray(self.func_vector + self.failed_evaluations))
        logger.info('Surrogate with failed points created! ({} points) '.format(self.curr_points))

    def _is_near_to_population(self, x: np.array, population: List[np.array] = None) -> bool:
        """
        Checks whether the next point is "near" to the already existing points
        Returns True if it is near to already existing points

        :param x: location to check nearness
        :param population : population to check nearness from, by default it will be set to current total population
        :return: bool: True/False
        """
        if population is None:
            population = self.total_population

        for member in population:
            dist = ut.distance(member, x)
            # estimating that if the distance is less that 1% of length between upper and lower bounds,
            # then the rows will be very close
            # TODO: THIS CAN LEAD TO PROBLEMS IF THETA BECOMES TOO LARGE
            if dist < 1e-2 * np.linalg.norm(self.l_bound-self.u_bound):
                return True
        return False

    def set_new_kernel_parameters(self, theta: Union[float, np.array] = 1.0, theta0: float = 1.0) -> None:
        """
        sets the kernel parameters: called from the user script.
        The covariance matrix needs to be re-calculated when this function is called

        :param theta: scaling of dimensions used with the distance function
        :param theta0: scaling of kernel function
        :return: None
        """
        pre_existing_points = self.curr_points
        logger.debug('Population size before kernel parameters update: {}'.format(self.curr_points))
        self.kernel.theta0 = theta0
        self.kernel.theta = theta
        # recalculate covariance matrix using add_point
        curr_pop = self.total_population.copy()
        curr_failed = self.failed_query_points.copy()
        curr_y = self.func_vector.copy()
        self.total_population = []
        self.failed_query_points = []
        self.func_vector = []
        self.curr_points = 0
        self.K = np.identity(self.curr_points)
        for pop, func in zip(curr_pop, curr_y):
            self.add_point(pop, func, if_check_nearness=True, is_failed=False)
        for pop in curr_failed:
            self.add_point(pop, if_check_nearness=True, is_failed=True)
        self.update_surrogate()
        logger.debug('Kernel function parameters changed: theta: {}, theta0: {}'.format(self.kernel.theta,
                                                                                        self.kernel.theta0))
        logger.debug('Population size after kernel parameters update: {}'.format(self.curr_points))
        logger.debug('Covariance matrix recalculated!')
        assert pre_existing_points == self.curr_points, 'Not all points are added to the krig surface '\
                                                        'after setting new kernel parameters'

    # getter functions for user access
    def get_total_population(self) -> Tuple[List[np.array], List[float]]:
        """
        returns the total population : for visualization purposes

        :return: total_population and function values at that location
        """
        return self.total_population, self.func_vector

    def get_current_iteration(self) -> int:
        """
        Returns the current iteration number

        :return: integer
        """
        return self.curr_iter

    def get_current_best(self) -> tuple:
        """
        returns the location of current best and value of current best

        :return: location, func_value
        """
        best_index = np.argmin(self.func_vector)
        return self.total_population[best_index], self.func_vector[best_index]

    def _require_all_cost_function(self, xs: List[np.array]) -> List[Tuple[np.array, float]]:
        """
        Function that evaluates all cost function, instead of a required fraction

        :param xs: list of locations to evaluate
        :return: list of tuples of locations and cost functions
        """
        pending = xs
        completed = list()
        failed = list()
        while len(pending) > 0:
            new_completed, pending, new_failed = self.cost_function(pending)
            completed += new_completed
            failed += new_failed

        if len(failed) > 0:
            warnings.warn('Could not evaluate all required points (' + str(len(failed)) + ' evaluations failed).',
                          category=RuntimeWarning)
        return completed

    def export_csv(self, path: str, **kwargs) -> None:
        """
        Writes the data from all completed function evaluations (i.e. x, y pairs) to a CSV file.
        The file starts with the following header row: x0, x1, ..., x[n_dim-1], y.
        Additional arguments will be forwarded to the csv.writer function, which can be used to control the
        formatting of the CSV file. To write a TSV instead of CSV, you can pass dialect='excel-tab'.

        :param path: Path to write the file to. If the file already exists, it will be overwritten.
        :param kwargs: Forwarded to csv.writer(), can be used to override the CSV dialect and set formatting options.
        """
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f, **kwargs)

            # write header (x0, x1, ..., xN-1, y)
            writer.writerow(['x' + str(n) for n in range(self.n_dim)] + ['y'])

            # write data
            for x, y in zip(self.total_population, self.func_vector):
                writer.writerow(list(x) + [y])
