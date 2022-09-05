.. _example-4:

Kriging
=======

Kriging or Gaussian process regression is a method of interpolation for which the interpolated values are modeled by a
Gaussian process governed by prior covariances. In this example, we show how PARyOpt can be used to generate response
surfaces using available data.
As with the previous examples, we shall use the standard parabola as the underlying function to be approximated. There
are several ways to use PARyOpt for Kriging, one of which is shown here. This is possibly the easiest and cleanest way
to perform Kriging using PARyOpt.

**Data generation**: Since the underlying function is known, we shall generate data by invoking this function
at some random locations within the bounds and storing them in an external file. This is achieved through the
following snippet:


.. code-block:: python

   def create_data_csv(function: Callable, filename: str, l_bound: np.array, u_bound: np.array) -> None:
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

**Data assimilation** : PARyOpt provides an :py:meth:`add_point` to add external data manually. The user has to supply
the `x` location and any available `y` values to this method to add data to the :py:class:`BayesOpt` instance.
An example usage of :py:meth:`add_point` using the above generated data can be:

.. code-block:: python

   def load_from_csv(b_opt: BayesOpt, filename: str) -> BayesOpt:
    """
    load data from csv file and add to PARyOpt
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

Note that the user has to manually invoke :py:meth:`update_surrogate`. This is currently for efficiency purposes and hope
to be replaced in the upcoming versions.

Finally, since the user wants to add data manually and does not want the standard initialization required for bayesian
optimization, we provide a switch :py:data:`do_init` to turn off the initialization. Since there is no cost function to
be optimized, the evaluator should be passed in an empty function for evaluation.

.. code-block:: python

    # dummy evaluators:
    evaluator = FunctionEvaluator(lambda x: 0.0)

    krig = BayesOpt(cost_function=evaluator,
                    l_bound=l_bound, u_bound=u_bound, n_dim=1,
                    n_init=0, do_init=False,        # ensures that initialization is not done.
                    kern_function='sqr_exp',
                    acq_func='LCB',
                    kappa_strategy=lambda curr_iter: 1000,
                    if_restart=False)
    krig = load_from_csv(krig, data_filename)
    krig.estimate_best_kernel_parameters(theta_bounds=[[0.001, 10.0]])

Note that since we are not providing any actual cost function here, :py:meth:`update_iter` does nothing useful. In case
the user is looking for an instantaneous Kriging model, i.e., creating a Kriging surface and updating it, the actual
cost function should be provided. Just like the previous examples, one may use :py:class:`evaluator` from example-3_
and do Kriging similar to optimization.

Now that the surrogate is created and hyper-parameters optimized, one can start querying it using
:py:meth:`evaluate_surrogate_at`

.. code-block:: python

   location = np.array([1.0])
   mean, variance = krig.evaluate_surrogate_at(location)



