.. _example-2:

Restart from previous state
===========================

This example shows how the user can deal with optimization failures, either due to hardware failure or due to wrong
selection of accessory functions. This will also be useful in cases of changing optimization platform but resuming the
same optimization task.

In order to create an optimization state, we shall first run a standard optimization problem for a certain number of
iterations. This can be done either by running example 0(:ref:`example-0`), or by re-doing the whole procedure. For the benefit
of the user, we shall take the latter way.

Restart
-------

In example 1, we have seen that our custom initialization is not the best compared to the default latin hypercube
sampling. In fact, such an example provides the best motivation for a restart. Hence, in this example, we provide
a custom initialization method, the same as in example 1 (:ref:`example-1`). Once the optimization is done for 10
iterations, we shall create another instance of :py:class:`BayesOpt` that starts from this existing optimization state.

.. code-block:: python

    restarted_bo = BayesOpt(cost_function=evaluator,
                            l_bound=l_bound, u_bound=u_bound, n_dim=n_dim,
                            n_init=2,
                            kern_function='sqr_exp',
                            acq_func='LCB',
                            kappa_strategy=my_kappa,
                            if_restart=True, restart_filename='opt_state.dat')

Note that the restarted optimization need not have the same accessory functions, like kernel, acquisition and kappa
strategy. By enabling ``if_restart`` and providing the ``restart_filename``, the framework re-creates the optimization
state from which the user can continue the optimization, for example,

.. code-block:: python

   restarted_bo.update_iter(5)

will update 5 iterations at once. This API helps to reduce redundant loops in the user code.

Intra-iteration restart
-----------------------

It has to be noted here that the evaluator also has an inbuilt check-pointing per iteration, so that hard interrupts
such as the ``KeyboardInterrupt`` can also be handled for restart. The user need not do any extra changes to enable
this **intra-iteration** restart functionality.

