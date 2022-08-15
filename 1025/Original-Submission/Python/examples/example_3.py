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
Example 3: Async Local evaluator
2 - dimensional

What you will learn:
* multiple optima per iteration
* Instantiating and using asynchronous evaluator
* extension to remote evaluators should be read from the docs
"""
from typing import List, Any
import sys
import os
import logging
import time
import numpy as np

from PARyOpt import BayesOpt
from PARyOpt.evaluators import FunctionEvaluator, AsyncLocalEvaluator, AsyncSbatchEvaluator
import examples.examples_all_functions as exf


# for AsyncLocalFunctionEvaluator
def folder_generator(directory, x) -> None:
    """
    prepares a given folder for performing the simulations. The cost function (out-of-script) will be executed
    in this directory for location x. Typically this involves writing a config file, generating/copying meshes and

    In our example, we are running a simple case and so does not require any files to be filled. We shall pass the
    location of cost function as a command line argument
    :param directory:
    :param x:
    :return:
    """
    with open(os.path.join(directory, 'config.txt'), 'w') as f:
        pass  # write file
    pass


def run_cmd(directory, x) -> List[Any]:
    """
    Command to run on local machine to get the value of cost function at x, in directory.
    In this example, we shall run the script example3_evaluator.py with the location as an argument.
    :param directory:
    :param x:
    :return:
    """
    eval_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'example3_evaluator.py')
    return [sys.executable, eval_path] + list(x[:])


def result_parser(directory, x) -> float:
    """
    Parses the result from a file and returns the cost function.
    The file is written be the actual cost function. One can also do post processing in this function and return the
    subsequent value. Based on the construct of our cost function example3_evaluator.py, the generated result.txt
    will be in this 'directory'
    :param directory:
    :param x:
    :return:
    """
    with open(os.path.join(directory, 'result.txt'), 'r') as f:
        return float(f.readline())


def user_defined_kappa(curr_iter, freq, t_const):
    """
    user defined kappa for multiple acquisition functions
    :param curr_iter:
    :param freq:
    :param t_const:
    :return:
    """
    kappa = 40.5 * (np.sin(curr_iter * np.pi / freq) + 1.5) * np.exp(-t_const*curr_iter)
    return kappa


if __name__ == "__main__":

    # logging setup
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler('../logs/example3_log_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S")), mode='a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.info('Performing a simple bayesian optimization on parabolic cost function with local asynchronous evaluator')

    # define basic parameters
    # lower and upper bounds
    l_bound = np.asarray([-4.0, -4.0])
    u_bound = np.asarray([10.0, 10.0])
    # dimensionality
    n_dim = 2
    # optima per iteration
    n_opt = 2
    # max number of iterations
    iter_max = 8

    jobs_dir = os.path.join(os.getcwd(), 'temp/opt_jobs')
    # parallel, asynchronous, out-of-script
    evaluator = AsyncLocalEvaluator(job_generator=folder_generator,
                                    run_cmd_generator=run_cmd,
                                    parse_result=result_parser,
                                    required_fraction=0.5, jobs_dir=jobs_dir)

    logger.info('Optimization evaluations are done in {} directory'.format(jobs_dir))
    # generate a list of kappa strategies (functions) that correspond to each acquisition function
    my_kappa_funcs = []
    for j in range(n_opt):
        my_kappa_funcs.append(lambda curr_iter_num, freq=10.*(j*j+2), t_const=0.8/(1. + j):
                              user_defined_kappa(curr_iter_num, freq=freq, t_const=t_const))

    b_opt = BayesOpt(cost_function=evaluator,
                     n_dim=n_dim, n_opt=n_opt, n_init=2,
                     u_bound=u_bound, l_bound=l_bound,
                     kern_function='matern_52',
                     acq_func='LCB', kappa_strategy=my_kappa_funcs,
                     if_restart=False)
    logger.info('BO initialized')

    for curr_iter in range(iter_max):
        b_opt.update_iter()
        if not curr_iter % 2:
            b_opt.estimate_best_kernel_parameters(theta_bounds=[[0.01, 10]])
        exf.visualize_fit(b_opt)

    # export cost function evaluations to a CSV file
    b_opt.export_csv(os.path.join(os.getcwd(), 'temp', 'data.csv'))

    # visualization of constructed surrogate
    exf.visualize_fit(b_opt)
    # get current best evaluated value
    best_location, best_value = b_opt.get_current_best()

    result_txt = 'Optimization done for {} iterations, best evaluation is at {} with cost: {}'. \
        format(b_opt.get_current_iteration(), best_location, best_value)

    logger.info(result_txt)
    print(result_txt)
    logger.info('Asynchronous bayesian optimization completed!')
