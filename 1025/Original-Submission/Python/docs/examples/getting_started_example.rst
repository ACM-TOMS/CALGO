.. _example-0:

Getting Started
===============

This example shows how one gets started with the optimization software.
One is expected to see through **example_0.py** for a better understanding.


Just like any optimization, we should have information about the following:

1. dimensionality of the problem -- ``n_dim``
2. cost function that has to be **minimized** -- ``function``
3. bounds on the variables   -- ``l_bound``, ``u_bound``

Some more parameters for Bayesian Optimization

4. number of initial evaluations for constructing the prior --``n_init``
5. type of kernel function -- ``kern_function``
6. number of evaluations per iteration -- ``n_opt``
7. platform of evaluation - local computer / remote computer
8. parallel / serial evaluations -- ``cost_function``
9. asynchronocity of evaluations -- ``cost_function``
10. acquisition function (list) -- ``acq_func``

    * parallelization of acquisition function -- ``kappa_strategy``

The ``kappa_strategy`` defines exploration vs exploitation of the optimizer. Those with a *large* ``kappa`` value
will explore and a *small* ``kappa`` value will exploit.

In this example, we shall solve a simple parabolic cost function, on a local machine with no parallelization. The evaluations
will, therefore, be fully synchronous. As part of the example, we shall do an exploration dominated search in the
optimization. The parameters that will be used are as follows:

.. code-block:: python

    from PARyOpt.evaluators import FunctionEvaluator
    import numpy as np
    from PARyOpt import BayesOpt

    n_dim = 1
    l_bound = np.asarray([-12.])
    u_bound = np.asarray([12.])

    n_init = 2
    kern_function = 'sqr_exp'       # squared exponential
    acq_func = 'LCB'                # lower confidence bound
    def my_cost_function(x: np.array) -> float:
            y = np.sum((x-2.5) ** 2 + 5)
            return float(y)
    # instantiate an evaluator that evaluates serially on the local machine
    evaluator = FunctionEvaluator(my_cost_function)
    def my_kappa(curr_iter: int) -> float:
        return 1000.0           # large value for exploration


Initialization
--------------

Having defined these parameters, we shall now initialize the optimizer:

.. code-block:: python

    b_opt = BayesOpt(cost_function=evaluator,
                     l_bound=l_bound, u_bound=u_bound, n_dim=n_dim,
                     n_init=2,
                     kern_function='sqr_exp',
                     acq_func='LCB',
                     n_opt=1,                           # default setting
                     kappa_strategy=my_kappa,
                     if_restart=False)

The stage is now set for optimization to be performed. Since this package does not provide any standard termination criteria,
the user is expected to design a termination based on the nature of the problem. In this example, we shall look at
a very simple termination criterion of number of iterations.

.. code-block:: python

   max_iter = 10

Update
------

The user shall manually update the optimization every iteration. This provides ways to post-process user required metrics
every iteration, as well as do a regular *hyper-parameter optimization* for optimized surrogate.

.. code-block:: python

    for curr_iter in range(max_iter):
        b_opt.update_iter()

Hyper parameter optimization
----------------------------

An implementation of the standard hyper parameter optimization is done in :py:meth:`estimate_best_kernel_parameters()`.
This minimizes a maximum likelihood estimate of the constructed surrogate and eventually sets the *optimal* kernel
length scale. It can be invoked by calling:

.. code-block:: python

   theta_min = 0.01
   theta_max = 50.
   b_opt.estimate_best_kernel_parameters(theta_bounds=[[theta_min, theta_max]])

Hyper parameter optimization need not be performed every iteration as the surrogate may not change much with the
addition of a single data point. Hence its call can be periodic based on the iteration number.


Surrogate query
---------------

Having constructed the surrogate, one may need to query it for several purposes, including visualization, post-processing
and termination criteria. This functionality is provided through the :py:meth:`evaluate_surrogate_at()` function. It
returns the value of the mean and variance of the surrogate at the queried location.

.. code-block:: python

    location_to_query = np.asarray([0.5])
    mean, variance = b_opt.evaluate_surrogate_at(location_to_query)


Logging
-------

PARyOpt uses the python `logging <https://docs.python.org/3/library/logging.html>`_ module for logging. The user has to
instantiate the logger in the main code. If the logger is not initiated, the logs will be streamed to ``stdout``. An
example of using the logger is also in the above example:

.. code-block:: python

   import logging, time

   logger = logging.getLogger()
   logger.setLevel(logging.INFO)  # either NOTSET, INFO, DEBUG, WARNING, ERROR, CRITICAL -- different levels of log
   log_file_name = 'example0_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S"))
   fh = logging.FileHandler(log_file_name, mode='a')
   # logging format
   formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
   fh.setFormatter(formatter)
   logger.addHandler(fh)

Saving data
-----------

The framework provides multiple ways to save data, particularly with ready methods to export in ``.csv`` format. It can be
done be calling:

.. code-block:: python

   b_opt.export_csv('my_data.csv')


Alternately, data can be custom exported, as ``get`` methods exist to get the population and the respective function values.

.. code-block:: python

   total_population, function_values = b_opt.get_total_population()


The next example shows how to custom change the various functions used in the optimization method.