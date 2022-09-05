.. _example-1:

Custom functions for surrogate construction
===========================================

This example discusses the modularity provided by the framework to change various optimizer and surrogate related
methods. More specifically, we look at changing the following accessory functions:

* :ref:`ex1-init`
* :ref:`ex1-kernel`
* :ref:`ex1-acq`
* :ref:`ex1-acq-opt`

Through this tutorial, we shall not solve any new problem. However, explanations will be provided to how can one add
custom definitions to these standard methods. We shall use an external implementation of the standard methods for the
user to easily compare with :ref:`example-0`. While the example shows the usage of all of these at once, the user is
encouraged to understand the effect of changing each of these accessory functions separately:

.. _ex1-init:

Initialization strategy
-----------------------

.. external constraints
.. known important locations
.. different sampling strategy

Several situations arise when the user has to specify his own initialization strategies. Some of these include:

* **Constraints**, such as equality, inequality and PDE based constraints. In such situations, all locations determined
  by the upper and lower bounds may not be feasible. PARyOpt requires all the initial evaluation points to be successful
  evaluations, hence a general latin hypercube sampling may not work in all cases.

* **Biased sampling** of search domain. There could be situations where the user knows strategic locations and thus can
  help PARyOpt to sample in those locations right from the beginning

* **User defined sampling strategy**. The user may want a completely different sampling strategy which has a better
  coverage of the domain for that specific optimization problem.

In order to specify custom initialization strategy, the user should define a function with the following signature:

.. code-block:: python

   def my_init_strategy(n_dim: int, n_init: int, l_bound: np.array, u_bound: np.array) -> List[np.array]

and then pass it during initialization to the parameter :code:`init_strategy=my_init_strategy`.


.. _ex1-kernel:

Kernel functions (Class)
------------------------

.. continuity and differentiability
.. region of influence prior information
.. periodic functions

Kernel function embeds most of the information about the **continuity** of the underlying function. Thus, it is one of the
critical parameters to be selected by the user. For example, Matern class of kernel functions have only finite
differentiability while the squared exponential kernel function is infinitely differentiable. Some other situations
where the user may want to change the kernel function are :

* **Underlying function continuity and differentiability**, as discussed above.
* **Prior information about periodicity** of the underlying function. In such situations, periodicity can be embedded
  into the kernel function. This can drastically improve the number of function evaluations. One such example of a
  periodic kernel function is:

.. math::
   k_{periodic}(x_1, x_2) = \theta_0 exp(-2 sin(\frac{\pi}{p} \frac{||x_1 - x_2||}{l})^2)

where :math:`p` is the periodicity interval and :math:`l` is the length scale of the kernel.

* **Anisotropic kernels** may be used when the function behaviour is different across different optimization parameters.
  For example, if the function varies logarithmically in one parameter and linearly in another parameter

To specify custom kernel functions, the user has to derive the base kernel class. An example to specify the standard
Matern 3/2 function by the user is as follows:

.. code-block:: python

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

This derived class requires defining the :code:`eval()` and :code:`derivative()` methods of the class. The parameter
:code:`theta` should contain all the hyper-parameters, for ex. length scale, related to the supplied kernel. This will
be used by PARyOpt during hyper-parameter optimization. This kernel class can be passed to the constructor through
the parameter :code:`kern_function=MyKernelFunction()`.


.. _ex1-acq:

Acquisition functions
---------------------

.. generally not required to change
.. safety constraints and other information about 'mean' value...

Acquisition function define informative regions of the surrogate. Since the software is designed for minimization, areas
of high information should have small acquisition values. In normal circumstances, the user should not need to change
this. However, if the user wants implementing **safety constraints** and **biased informativeness**, it can be be done
by passing a user defined function, with the following signature, to the constructor as
:code:`acq_func=my_acquisition_function`.

.. code-block:: python

   def my_acquisition_function(mean: float, variance: float, curr_best: float = 0., kappa: float = 1.) -> float:

.. _ex1-acq-opt:

Acquisition optimizer
---------------------

.. dimensionality
.. in-house optimizer
.. GA/PSO type heuristic optimizers

Bayesian optimization proceeds by evaluating locations with maximum information content. Hence, it requires finding the
optimum (minimum in PARyOpt) of the acquisition function (cheap to evaluate). While the core software comes with a
standard ``Powell`` algorithm, the user may quite often want to change this for a better global optimum. Some of the
reasons include:

* **Dimensionality of optimization** could impose restrictions on the type of optimizer used
* **Heuristic optimization** can be an alternative for multi-modal functions. Since evaluating the acquisition is fast,
  these optimizers will also be efficient.
* Robust **in-house optimizers** may be available with research groups tailored for specific problems.

The process of adding user-defined acquisition optimizer is very similar to defining custom acquisition function. The
function signature is :

.. code-block:: python

   def my_acq_optimizer(func: Callable[[np.array], float], x0: np.array, l_bound: np.array, u_bound: np.array) -> np.array:

which is passed to the constructor as :code:`acq_func_optimizer=my_acq_optimizer`.
