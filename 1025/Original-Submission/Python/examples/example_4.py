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
Example 4: Kriging using
1 - dimensional

What you will learn:
* saving and retrieving data from csv files into bayes opt
* kriging
"""
import csv
import logging
import time

import numpy as np
from typing import Callable

from PARyOpt import BayesOpt
import examples.examples_all_functions as exf
from PARyOpt.evaluators import FunctionEvaluator


def load_from_csv(b_opt: BayesOpt, filename: str) -> BayesOpt:
    """
    load data from csv file and add to PARyOpt
    :param b_opt:
    :param filename:
    :return:
    """
    with open(filename, 'r') as csvfile:
        csv_file_lines = csv.reader(csvfile, delimiter=',')
        for row_num, row in enumerate(csv_file_lines):
            if row_num == 0:
                # skipping the header
                pass
            else:
                b_opt.add_point(x=np.asarray([float(row[0])]), y=float(row[-1]),
                                if_check_nearness=True)
    b_opt.update_surrogate()

    return b_opt


def create_data_csv(function: Callable, filename: str, l_bound: np.array, u_bound: np.array) -> None:
    """
    take the cost function and save data to a csv file
    :param l_bound:
    :param u_bound:
    :param function:
    :param filename:
    :return:
    """
    # generate some random locations -- 7
    normalized_population = np.random.ranf((7, ))
    real_population = l_bound + normalized_population * (u_bound - l_bound)
    real_population = [np.asarray([p]) for p in real_population]
    # evaluate the values
    real_functions = [float(function(p)) for p in real_population]

    # write into file
    with open(filename, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow('x y')
        for x, y in zip(real_population, real_functions):
            writer.writerow(list(x) + [y])


if __name__ == "__main__":

    # set up logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler('../logs/example4_log_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S")), mode='a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.info('Performing a simple bayesian optimization on parabolic cost function with local asynchronous evaluator')

    # bounds
    l_bound = np.asarray([-12.])
    u_bound = np.asarray([12.])
    # dummy evaluators:
    evaluator = FunctionEvaluator(lambda x: 0.0)
    # set up a BayesOpt object for kriging
    krig = BayesOpt(cost_function=evaluator,
                    l_bound=l_bound, u_bound=u_bound, n_dim=1,
                    n_init=0, do_init=False,        # ensures that initialization is not done.
                    kern_function='sqr_exp',
                    acq_func='LCB',
                    kappa_strategy=lambda curr_iter: 0.01,
                    if_restart=False)

    # Stage 1: get some data in a csv file. Usually the user has data generated from an external source,
    # but in this example, we shall generate and store the data
    # create data in a csv file
    data_filename = '../temp/example_4_data.csv'
    create_data_csv(exf.parabolic_cost_function, data_filename, l_bound, u_bound)

    # Stage 2: add points into krig from the file
    krig = load_from_csv(krig, data_filename)

    # first optimize the hyperparameters
    krig.estimate_best_kernel_parameters(theta_bounds=[[0.001, 10.0]])

    # surrogate is now ready for query -- lets visualize it
    exf.visualize_fit(krig)

    logger.info('Kriging surrogate is ready!')
