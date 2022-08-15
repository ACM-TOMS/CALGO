Examples
========

This provides a demonstration of the exhaustive functionality of PARyOpt. These are very intricately connected to
the examples on BitBucket, so the user is suggested to go through them simultaneously. The examples are structured as
follows:


.. csv-table:: Description of Examples
   :header: "Example #", "Description"
   :widths: 15, 30
   :align: right

   "Example **0**

   :ref:`example-0`", "A dive into setting up the optimization routine."
   "Example **1**

   :ref:`example-1`", "Using the same problem as example 0, the modularity of the
   framework is demonstrated. All the functions that can be customized are shown."
   "Example **2**

   :ref:`example-2`", "Explains the restart capabilities of the framework and how
   one can use in case of (hardware/resource) failure."
   "Example **3**

   :ref:`example-3`", "Asynchronocity is introduced and a local asynchronous implementation
   of example 0 is shown. This implementation can be easily extended to an asynchronous
   remote evaluator."
   "Example **4**

   :ref:`example-4`", "Explains the kriging functionality of the framework and how data
   can be assimilated to perform kriging"



.. toctree::
   :maxdepth: 3
   :caption: Tutorials

   getting_started_example
   user_defined_function_example
   restart_example
   async_local_evaluator_example
   kriging_example
