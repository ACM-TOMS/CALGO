.. _sgtelib:

Surrogate Library
========================

The *SGTELIB* library is a dynamic surrogate modelling library. It is used in the *Search* step of Mads to dynamically construct models from the previous evaluations.
During a *Search* step that uses *SGTELIB*, models of the objective and the constraints are constructed and a surrogate subproblem involving these models is optimized.
The resulting solutions are the next candidates for evaluation by the true problem.

| Models from the *SGTELIB* library can be used by setting the parameter ``SGTELIB_MODEL_SEARCH`` to ``yes`` or ``true``.


Models
-------------------

Models in *SGTELIB* are defined by using a succession of field names and field values.
To choose a model, the parameter ``SGTELIB_MODEL_DEFINITION`` must be used followed by the field name ``TYPE``, and then by the model type.
The subsequent fields enable to define the settings of the model.
Each field name is made of one single word and each field value is made of one single word or numerical value.

Example : ``SGTELIB_MODEL_DEFINITION TYPE <model type> FIELD1 <field 1 value> FIELD2 <field 2 value>``

The section below describes the models and settings available.


Types of models
""""""""""""""""""""""

Below is the list of all possible models and their authorized fields.

.. _prs:

``PRS``
""""""""
| PRS (Polynomial Response Surface) is a type of model.
| Authorized fields:

* :ref:`degree` (Can be optimized)

* :ref:`ridge` (Can be optimized)

* :ref:`budget`: Defines the budget allocated for parameter optimization.

* :ref:`output`: Defines the output text file.

| Examples:
| ``TYPE PRS DEGREE 2``
| ``TYPE PRS DEGREE OPTIM RIDGE OPTIM``


.. _prs_edge:

``PRS_EDGE``
""""""""""""""

| PRS_EDGE (Polynomial Response Surface EDGE) is a type of model that allows to model discontinuities at :math:`0` by using additional basis functions.
| Authorized fields:

* :ref:`degree` (Can be optimized)

* :ref:`ridge` (Can be optimized)

* :ref:`budget`: Defines the budget allocated for parameter optimization.

* :ref:`output`: Defines the output text file.

| Examples:
| ``TYPE PRS_EDGE DEGREE 2``
| ``TYPE PRS_EDGE DEGREE OPTIM RIDGE OPTIM``


.. _prs_cat:

``PRS_CAT``
""""""""""""""
| PRS_CAT (Categorical Polynomial Response Surface) is a type of model that allows to build one PRS model for each different value of the first component of :math:`x`.
| Authorized fields:

* :ref:`degree` (Can be optimized)
* :ref:`ridge` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE PRS_CAT DEGREE 2``
| ``TYPE PRS_CAT DEGREE OPTIM RIDGE OPTIM``


.. _rbf:

``RBF``
""""""""""""""
| RBF (Radial Basis Function) is a type of model.
| Authorized fields:

* :ref:`kernel_type` (Can be optimized)
* :ref:`kernel_shape` (Can be optimized)
* :ref:`distance_type` (Can be optimized)
* :ref:`ridge` (Can be optimized)
* :ref:`preset`: Defines the type of RBF model used.
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE RBF KERNEL_TYPE D1 KERNEL_SHAPE OPTIM DISTANCE TYPE NORM2``


.. _ks:

``KS``
""""""""""""""
| KS (Kernel Smoothing) is a type of model.
| Authorized fields:

* :ref:`kernel_type` (Can be optimized)
* :ref:`kernel_shape` (Can be optimized)
* :ref:`distance_type` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM`` 


.. _kriging:

``KRIGING``
""""""""""""""
| KRIGING is a type of model.
| Authorized fields:

* :ref:`ridge` (Can be optimized)
* :ref:`distance_type` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE KRIGING``


.. _lowess:

``LOWESS``
""""""""""""""
| LOWESS (Locally Weighted Regression) is a type of model (from [TaAuKoLed2016]_).
| Authorized fields:

* :ref:`degree`: Must be 1 (default) or 2 (Can be optimized).
* :ref:`ridge` (Can be optimized)
* :ref:`kernel_type` (Can be optimized)
* :ref:`kernel_shape` (Can be optimized)
* :ref:`distance_type` (Can be optimized)
* :ref:`preset`: Defines how the weight of each data point is computed.
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE LOWESS DEGREE 1``
| ``TYPE LOWESS DEGREE OPTIM KERNEL_SHAPE OPTIM KERNEL_TYPE D1``
| ``TYPE LOWESS DEGREE OPTIM KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM DISTANCE TYPE OPTIM``


.. _cn:

``CN``
""""""""""""""
| CN (Closest Neighbours) is a type of model.
| Authorized fields:

* :ref:`distance_type` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE CN``


.. _ensemble:

``ENSEMBLE``
""""""""""""""
| ENSEMBLE is a type of model that uses multiple models simultaneously.
| Authorized fields:

* :ref:`weight`: Defines how the ensemble weights are computed.
* :ref:`metric`: Defines which metric is used to compute the weights.
* :ref:`distance_type`: This parameter is transfered to the models contained in the Ensemble.
* :ref:`preset`: Defines the selection of models in the ensemble.
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE ENSEMBLE WEIGHT SELECT METRIC OECV``
| ``TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV DISTANCE TYPE NORM2 BUDGET 100``


.. _ensemble_stat:

``ENSEMBLE_STAT``
""""""""""""""""""

| ENSEMBLE_STAT is a type of model (from [AuLedSa2021]_).
| Authorized fields:

* all the fields from :ref:`ensemble` (with different default values though).
* :ref:`uncertainty`: Selects an alternative for the uncertainty (smooth or nonsmooth).
* :ref:`size_param`: Defines the size parameter (different meaning depending on the value of UNCERTAINTY).
* :ref:`sigma_mult`: Defines the scaling factor of the uncertainty.
* :ref:`lambda_p`: Defines the shape parameter of the probability of feasibility.
* :ref:`lambda_pi`: Defines the shape parameter of the probability of improvement.

| Example:
| ``TYPE ENSEMBLE_STAT UNCERTAINTY SMOOTH WEIGHT SELECT5 METRIC RMSECV SIZE_PARAM 15``



The following table summarizes the possible fields for every model.

.. csv-table:: Model authorized fields
   :header: "Model type", :ref:`degree`, :ref:`ridge`, :ref:`kernel_type`, :ref:`kernel_shape`, :ref:`distance_type`, :ref:`preset`, :ref:`weight`, :ref:`metric`, :ref:`uncertainty`,:ref:`budget`, :ref:`output`

   :ref:`prs`,          ✔,  ✔,  ,    ,    ,   ,  ,  ,  , ✔, ✔
   :ref:`prs_edge`,     ✔,  ✔,  ,    ,    ,   ,  ,  ,  , ✔, ✔
   :ref:`prs_cat`,      ✔,  ✔,  ,    ,    ,   ,  ,  ,  , ✔, ✔
   :ref:`rbf`,           ,  ✔,  ✔,  ✔,  ✔, ✔,   ,  ,  , ✔, ✔
   :ref:`ks`,            ,   ,  ✔,  ✔,  ✔,   ,   ,  ,  , ✔, ✔
   :ref:`kriging`,       ,  ✔,  ,    ,   ✔,  ,    ,  ,  , ✔, ✔
   :ref:`lowess`,       ✔, ✔, ✔,  ✔,   ✔, ✔,    ,  ,  , ✔, ✔
   :ref:`cn`,            ,   ,  ,    ,   ✔,  ,    ,  ,  , ✔, ✔
   :ref:`ensemble`,      ,   ,  ,    ,   ✔, ✔,  ✔, ✔,  , ✔, ✔
   :ref:`ensemble_stat`, ,   ,  ,    ,   ✔, ✔,  ✔, ✔, ✔, ✔, ✔


Main model parameters
""""""""""""""""""""""""""

Below is the list of fields and their descriptions.

.. _degree:

``DEGREE``
""""""""""""""
| The field name DEGREE defines the degree of a polynomial response surface. The value must be an integer :math:`\geq 1`.
| Allowed for models of type: :ref:`prs`, :ref:`prs_edge`, :ref:`prs_cat` and :ref:`lowess`.
| Default value: 5

* For PRS models, the default degree is 2.
* For LOWESS models, the degree must be 1 (default) or 2.

| Example:
| ``TYPE PRS DEGREE 3 defines a PRS model of degree 3.``
| ``TYPE PRS_EDGE DEGREE 2 defines a PRS_EDGE model of degree 2.``
| ``TYPE LOWESS DEGREE OPTIM defines a LOWESS model where the degree is optimized.``


.. _ridge:

``RIDGE``
""""""""""""""
| The field name RIDGE defines the regularization parameter of the model.
| Allowed for models of type: :ref:`prs`, :ref:`prs_edge`, :ref:`prs_cat`, :ref:`rbf`, :ref:`kriging` and :ref:`lowess`.
| Possible values: Real value :math:`\geq 0`. Recommended values are :math:`0` and :math:`0.001`.
| Default value: :math:`0.001`.

| Example:
| ``TYPE PRS DEGREE 3 RIDGE 0`` defines a PRS model of degree 3 with no ridge.
| ``TYPE PRS DEGREE OPTIM RIDGE OPTIM`` defines a PRS model where the degree and ridge coefficient are optimized.


.. _kernel_type:

``KERNEL_TYPE``
""""""""""""""""
| The field name KERNEL_TYPE defines the type of kernel used in the model. The field name ``KERNEL`` is equivalent.
| Allowed for models of type: :ref:`rbf`, :ref:`lowess` and :ref:`ks`.
| Possible values:

* ``D1``: Gaussian kernel
* ``D2``: Inverse Quadratic Kernel
* ``D3``: Inverse Multiquadratic Kernel
* ``D4``: Bi-quadratic Kernel
* ``D5``: Tri-cubic Kernel
* ``D6``: Exponential Sqrt Kernel
* ``D7``: Epanechnikov Kernel
* ``I0``: Multiquadratic Kernel
* ``I1``: Polyharmonic splines, degree 1
* ``I2``: Polyharmonic splines, degree 2
* ``I3``: Polyharmonic splines, degree 3
* ``I4``: Polyharmonic splines, degree 4
* ``OPTIM``: The type of kernel is optimized

| Default value: ``D1``, except for RBF models where it is ``I2``.

| Example:
| ``TYPE KS KERNEL_TYPE D2`` defines a KS model with Inverse Quadratic Kernel.
| ``TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM`` defines a KS model with optimized kernel shape and type.


.. _kernel_shape:

``KERNEL_SHAPE``
""""""""""""""""""
| The field name KERNEL_SHAPE defines the shape coefficient of the kernel function. The field name ``KERNEL_COEF`` is equivalent. Note that this field name has no impact for kernel types ``I1``, ``I2``, ``I3`` and ``I4`` because these kernels do not include a shape parameter.
| Allowed for models of type: :ref:`rbf`, :ref:`ks` and :ref:`lowess`.
| Possible values: Real value :math:`\geq 0`. Recommended range is :math:`[0.1; 10]`. For KS and LOWESS model, small values lead to smoother models.
| Default value: By default, the kernel coefficient is optimized.

| Example:
| ``TYPE RBF KERNEL_TYPE D4 KERNEL_SHAPE 10`` defines a RBF model with an inverse bi-quadratic kernel of shape coefficient :math:`10`.
| ``TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM`` defines a KS model with optimized kernel shape and type.


.. _distance_type:

``DISTANCE_TYPE``
""""""""""""""""""
| The field name DISTANCE_TYPE defines the distance function used in the model.
| Allowed for models of type: :ref:`rbf`, :ref:`ks`, :ref:`kriging`, :ref:`lowess`, :ref:`cn`, :ref:`ensemble` and :ref:`ensemble_stat`.
| Possible values:

* ``NORM1``: Euclidian distance
* ``NORM2``: Distance based on norm :math:`1`
* ``NORMINF``: Distance based on norm :math:`1`
* ``NORM2_IS0``: Tailored distance for discontinuity in :math:`0`
* ``NORM2_CAT``: Tailored distance for categorical models

| Default value: ``NORM2``.

| Example:
| ``TYPE KS DISTANCE NORM2_IS0`` defines a KS model tailored for VAN optimization.


.. _preset:

``PRESET``
""""""""""""""
| The field name PRESET defines the type of model used when applicable.
| Allowed for models of type: :ref:`rbf`, :ref:`lowess`, :ref:`ensemble` and :ref:`ensemble_stat`.

* When applied to :ref:`rbf` models, PRESET defines the type of RBF.
      Possible values:

      * ``O``: RBF with linear terms and orthogonal constraints
      * ``R``: RBF with linear terms and regularization term
      * ``I``: RBF with incomplete set of basis functions (see [AuKoLedTa2016]_ for RBFI models)

      |
      | Default value: ``I``.

      | Example:
      | ``TYPE RBF PRESET O``

* When applied to :ref:`lowess` models [TaAuKoLed2016]_, PRESET defines how the weight :math:`w_i` of each data point :math:`x_i` is computed.
      Possible values:

      * ``D``: :math:`w_i=\phi(d_i)` where :math:`\phi` is the kernel of type and shape defined by the fields :ref:`kernel_type` and :ref:`kernel_shape`, respectively, and :math:`d_i` is the distance between the prediction point and the data point :math:`x_i`
      * ``DEN``: :math:`w_i=\phi(d_i/d_q)` where :math:`d_q` is the distance between the prediction point and the :math:`q^{th}` closest data point, and :math:`d_q` is computed with an empirical method
      * ``DGN``: :math:`w_i=\phi(d_i/d_q)` where :math:`d_q` is computed with the Gamma method
      * ``RE``: :math:`w_i=\phi(r_i)` where :math:`r_i` is the rank of :math:`x_i` in terms of distance to the prediction point, and :math:`r_i` is computed with empirical method
      * ``RG``: :math:`w_i=\phi(r_i)` where :math:`r_i` is computed with the Gamma method
      * ``REN``: same as ``RE`` but the ranks are normalized in :math:`[0,1]`
      * ``RGN``: same as ``RG`` but the ranks are normalized in :math:`[0,1]`

      |      
      | Default value: ``DGN``.

      | Example:
      | ``TYPE LOWESS PRESET RE``

* When applied to :ref:`ensemble` or :ref:`ensemble_stat` models, PRESET determines the selection of models in the ensemble.
      Possible values:

      * ``DEFAULT``: selection of 18 models of types :ref:`prs`, :ref:`ks`, :ref:`rbf` and :ref:`cn` with various settings
      * ``KS``: selection of 7 models of type :ref:`ks` with various kernel shapes
      * ``PRS``: selection of 7 models of type :ref:`prs` with various degrees
      * ``IS0``: selection of 30 models of type :ref:`prs_edge`, :ref:`ks`, :ref:`rbf` with various settings and DISTANCE_TYPE set to NOMR2_IS0
      * ``CAT``: selection of 30 models of type :ref:`prs_edge`, :ref:`ks`, :ref:`rbf` with various settings and DISTANCE_TYPE set to NOMR2_CAT
      * ``SUPER1``: selection of 4 models of types :ref:`prs`, :ref:`ks`, :ref:`rbf` and :ref:`lowess`
      * ``SMALL``: selection of 3 models of types :ref:`prs`, :ref:`ks` and :ref:`rbf`

      |
      | Default value: ``DEFAULT``.

      | Example:
      | ``TYPE ENSEMBLE PRESET SUPER1``


.. _weight:

``WEIGHT``
""""""""""""""
| The field name WEIGHT defines the method used to compute the weights :math:`\boldsymbol{w}` of the ensemble of models. The field name ``WEIGHT_TYPE`` is equivalent.
| Allowed for models of type: :ref:`ensemble` and :ref:`ensemble_stat`.
| Possible values:

* ``WTA1``: :math:`w_k \propto \mathcal{E}_{sum} - \mathcal{E}_k`
* ``WTA3``: :math:`w_k \propto (\mathcal{E}_k + \alpha\mathcal{E}_{mean})^{\beta}`
* ``SELECT``: :math:`w_k \propto 1` if :math:`\mathcal{E}_k = \mathcal{E}_{min}` (only the best model is selected)
* ``SELECTN``: :math:`w_k \propto \mathcal{E}_{sum}^N - \mathcal{E}_k` (for :math:`N=1,2,\dots,6`)
* ``OPTIM``: :math:`\boldsymbol{w}` minimizes :math:`\mathcal{E}(\boldsymbol{w})`

Where :math:`\mathcal{E}_k` is the error metric (defined by the field name :ref:`metric`) of the :math:`k^{th}` model in the ensemble,
:math:`\mathcal{E}_{sum}` is the cumulated error of all models,
:math:`\mathcal{E}_{min}` is the minimal error,
:math:`\mathcal{E}_{mean}` is the average error,
:math:`\alpha=0.05`, :math:`\beta=-1`,
and :math:`\mathcal{E}_{sum}^N` is the cumulated error metric of the :math:`N` best models.

| Default value: ``SELECT`` for :ref:`ensemble` models, ``SELECT3`` for :ref:`ensemble_stat` models with :ref:`uncertainty` set to ``SMOOTH``, and  ``SELECT4`` for :ref:`ensemble_stat` models with :ref:`uncertainty` set to ``NONSMOOTH``.

| Example:
| ``TYPE ENSEMBLE WEIGHT SELECT METRIC RMSECV`` defines an ensemble of models which selects the model that has the best RMSECV.
| ``TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV`` defines an ensemble of models where the weights :math:`\boldsymbol{w}` are computed to minimize the RMSECV of the model.
| ``TYPE ENSEMBLE WEIGHT SELECT3 METRIC OECV`` defines an ensemble of models which selects the 3 models that have the best OECV.


.. _uncertainty:

``UNCERTAINTY``
"""""""""""""""
(specific to :ref:`ensemble_stat` models)

| The field name UNCERTAINTY defines the type of uncertainty used in ENSEMBLE_STAT models. 
| Possible values:

* ``SMOOTH``: Smooth alternative of the uncertainty (default)
* ``NONSMOOTH``: Nonmooth alternative of the uncertainty

| Example:
| ``TYPE ENSEMBLE_STAT UNCERTAINTY NONSMOOTH``


.. _size_param:

``SIZE_PARAM``
""""""""""""""""
(advanced parameter specific to :ref:`ensemble_stat` models)

| The field name SIZE_PARAM defines the size of the directions of either :

- the simplex used to compute the simplex gradients of the models if the field :ref:`uncertainty` is set to ``SMOOTH``
- the positive spanning set used to compare models values if the field :ref:`uncertainty` is set to ``NONSMOOTH``

| Possible values: Real value :math:`\geq 0`. Recommended range is :math:`[0.001; 0.1]`.
| Default value: :math:`0.001` if the field UNCERTAINTY is set to ``SMOOTH``, :math:`0.005` if the field UNCERTAINTY is set to ``NONSMOOTH``.

| Example:
| ``TYPE ENSEMBLE_STAT UNCERTAINTY SMOOTH SIZE_PARAM 0.003``


.. _sigma_mult:

``SIGMA_MULT``
""""""""""""""""
(advanced parameter specific to :ref:`ensemble_stat` models)

| The field name SIGMA_MULT defines the scaling factor of the uncertain to be multiplied by the variance of already sampled function values.

| Possible values: Real value :math:`\geq 0`. Recommended range is :math:`[1; 100]`.
| Default value: :math:`10`.

| Example:
| ``TYPE ENSEMBLE_STAT UNCERTAINTY NONSMOOTH SIGMA_MULT 30``


.. _lambda_p:

``LAMBDA_P``
""""""""""""""""
(advanced parameter specific to :ref:`ensemble_stat` models)

| The field name LAMBDA_P defines the shape parameter of the *probability of feasibility* (P).

| Possible values: Real value :math:`\geq 0`. Recommended range is :math:`[0.1; 10]`.
| Default value: :math:`3` if the field UNCERTAINTY is set to ``SMOOTH``, :math:`1` if the field UNCERTAINTY is set to ``NONSMOOTH``.

| Example:
| ``TYPE ENSEMBLE_STAT UNCERTAINTY NONSMOOTH LAMBDA_P 1.5``


.. _lambda_pi:

``LAMBDA_PI``
""""""""""""""""
(advanced parameterspecific to :ref:`ensemble_stat` models)

| The field name LAMBDA_PI defines the shape parameter of the *probability of improvement* (PI).

| Possible values: Real value :math:`\geq 0`. Recommended range is :math:`[0.01; 3]`.
| Default value: :math:`0.1` if the field UNCERTAINTY is set to ``SMOOTH``, :math:`0.5` if the field UNCERTAINTY is set to ``NONSMOOTH``.

| Example:
| ``TYPE ENSEMBLE_STAT UNCERTAINTY NONSMOOTH LAMBDA_PI 0.3``


.. _output:

``OUTPUT``
""""""""""""""
Defines a text file in which model information are recorded. Allowed for ALL types of model.




Parameter optimization and selection
""""""""""""""""""""""""""""""""""""""""

Below is the list of some field names and values that influence the behaviour of other fields.

.. _optim:

``OPTIM``
""""""""""""""
| The field value OPTIM indicates that the model parameter must be optimized. The default optimization criteria is the AOECV error metric (except for ENSEMBLE_STAT models where it is OECV).
| Parameters that can be optimized:

* :ref:`degree`
* :ref:`ridge`
* :ref:`kernel_type`
* :ref:`kernel_shape`
* :ref:`distance_type`

| Example:
| ``TYPE PRS DEGREE OPTIM``
| ``TYPE LOWESS DEGREE OPTIM KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM METRIC ARMSECV``


.. _metric:

``METRIC``
""""""""""""""
| The field name METRIC defines the metric used to select the parameters of the model (including the weights of Ensemble models).
| Allowed for ALL types of model.
| Possible values:

* ``EMAX``: Error Max
* ``EMAXCV``: Error Max with Cross-Validation
* ``RMSE``: Root Mean Square Error
* ``RMSECV``: RMSE with Cross-Validation
* ``OE``: Order Error
* ``OECV``: Order Error with Cross-Validation [AuKoLedTa2016]_
* ``LINV``: Invert of the Likelihood
* ``AOE``: Aggregate Order Error
* ``AOECV``: Aggregate Order Error with Cross-Validation [TaAuKoLed2016]_

| Default value: ``AOECV``, except for :ref:`ensemble_stat` models where it is ``OECV``.

| Example:
| ``TYPE ENSEMBLE WEIGHT SELECT METRIC RMSECV`` defines an ensemble of models which selects the model that has the best RMSECV.


.. _budget:

``BUDGET``
""""""""""""""
| Budget for model parameter optimization. The number of sets of model parameters that are tested is equal to the optimization budget multiplied by the number of parameters to optimize.
| Allowed for ALL types of model.
| Default value: :math:`20`

| Example:
| ``TYPE LOWESS KERNEL_SHAPE OPTIM METRIC AOECV BUDGET 100``
| ``TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV BUDGET 50``




Surrogate subproblem formulations
-------------------------------------

The *SGTELIB* library offers different formulations of the surrogate subproblem to be optimized at the *Search* step (see [TaLeDKo2014]_).
The ``SGTELIB_MODEL_FORMULATION`` parameter enables to choose a formulation, and the parameter ``SGTELIB_MODEL_DIVERSIFICATION`` enables to adjust a diversification parameter.


``SGTELIB_MODEL_FORMULATION``
""""""""""""""""""""""""""""""

| The formulations of the surrogate subproblem involve various quantities.
| :math:`\hat f` denotes a model of the objective :math:`f` and :math:`\hat c_j` a model of the constraint :math:`c_j`, :math:`j=1,2,\dots,m`. For :math:`x\in X`, :math:`\sigma_f(x)` denotes the uncertainty associated with the prediction :math:`\hat f(x)`, and :math:`\sigma_j(x)` denotes the uncertainty associated with the prediction :math:`\hat c_j(x)`, :math:`j=1,2,\dots,m`. This uncertainty depends on the model chosen.

| For a :ref:`kriging` model, :math:`\sigma_f(x)` (or :math:`\sigma_j(x)`) is readily available through the standard deviation that the model natively produces.
| For an :ref:`ensemble_stat` model, the uncertainty is constructed by comparing the predictions of the ensemble models (see [AuLedSa2021]_).
| For any other model except ENSEMBLE, :math:`\sigma_f(x)` (or :math:`\sigma_j(x)`) is computed with the distance from :math:`x` to previously evaluated points.
| Finally, for an :ref:`ensemble` model, the uncertainty is computed through a weighted sum of the squared uncertainties of the ensemble models.

| There are eight different formulations that can be chosen with the parameter ``SGTELIB_MODEL_FORMULATION``. Some formulations involve a parameter :math:`\lambda` that is described later.

* ``FS`` (default):

.. math::

      \min_{x\in X}&\ \ \hat f(x)-\lambda\hat\sigma_f(x) \\
      \mathrm{s.t.}&\ \ \hat c_j(x)-\lambda\hat\sigma_j(x)\leq0,\ \ j=1,2,\dots,m

* ``FSP``:

.. math::

      \min_{x\in X}&\ \ \hat f(x)-\lambda\hat\sigma_f(x) \\
      \mathrm{s.t.}&\ \ \mathrm{P}(x)\geq 0.5

where :math:`\mathrm{P}` is the *probability of feasibility* which is the probability that a given point is feasible.

* ``EIS``:

.. math::

      \min_{x\in X}&\ -\mathrm{EI}(x)-\lambda\hat\sigma_f(x) \\
      \mathrm{s.t.}&\ \ \hat c_j(x)-\lambda\hat\sigma_j(x)\leq0,\ \ j=1,2,\dots,m

where :math:`\mathrm{EI}` is the *expected improvement* that takes into account the probability of improvement and
the expected amplitude thereof.

* ``EFI``:

.. math::
 
      \min_{x\in X}\ -\mathrm{EFI}(x)

where :math:`\mathrm{EFI}` is the *expected feasible improvement* : :math:`\mathrm{EFI} = \mathrm{EI}\times\mathrm{P}`.

* ``EFIS``:

.. math::
  
      \min_{x\in X}\ -\mathrm{EFI}(x)-\lambda\hat\sigma_f(x)

* ``EFIM``:

.. math::
  
      \min_{x\in X}\ -\mathrm{EFI}(x)-\lambda\hat\sigma_f(x)\mu(x)

where :math:`\mu` is the *uncertainty in the feasibility* : :math:`\mu = 4\mathrm{P}\times(1-\mathrm{P})`.

* ``EFIC``:

.. math::

      \min_{x\in X}\ -\mathrm{EFI}(x)-\lambda(\mathrm{EI}(x)\mu(x)
      +\mathrm{P}(x)\hat\sigma_f(x))

* ``PFI``:

.. math::
  
      \min_{x\in X}\ -\mathrm{PFI}(x)

where :math:`\mathrm{PFI}` is the *probability of improvement* : :math:`\mathrm{PFI} = \mathrm{PI}\times\mathrm{P}`,
with :math:`\mathrm{PI}` being the *probability of improvement* which is the probability that the objective decreases from the best known value at a given point.


| Example:
| ``SGTELIB_MODEL_DEFINITION TYPE KRIGING``
| ``SGTELIB_MODEL_FORMULATION EFIC``
| The two lines above define a surrogate subproblem based on the EFIC formulation that will involve kriging models.


``SGTELIB_MODEL_DIVERSIFICATION``
""""""""""""""""""""""""""""""""""

| The exploration parameter :math:`\lambda` enables to control the exploration of the search space against the intensification in the most promising areas. A higher :math:`\lambda` favors exploration whereas a lower :math:`\lambda` favors intensification.

| :math:`\lambda` is a real value in :math:`[0,1]` defined by the parameter ``SGTELIB_MODEL_DIVERSIFICATION``.
| Default value : :math:`0.01`.

| Example:
| ``SGTELIB_MODEL_DEFINITION TYPE ENSEMBLE``
| ``SGTELIB_MODEL_FORMULATION FSP``
| ``SGTELIB_MODEL_DIVERSIFICATION 0.1``
| The three lines above define a surrogate subproblem based on the FSP formulation with an exploration parameter equals to :math:`0.1` that will involve ensemble models.



.. topic:: References


  .. [TaAuKoLed2016] B.Talgorn, C.Audet, M.Kokkolaras and S.Le Digabel.
    Locally weighted regression models for surrogate-assisted design optimization.
    *Optimization and Engineering*, 19(1):213–238, 2018.
  
  .. [TaLeDKo2014] B.Talgorn, S.Le Digabel and M.Kokkolaras.
    Statistical Surrogate Formulations for Simulation-Based Design Optimization.
    *Journal of Mechanical Design*, 137(2):021405–1–021405–18, 2015
  
  .. [AuKoLedTa2016] C.Audet, M.Kokkolaras, S.Le Digabel and B.Talgorn.
    Order-based error for managing ensembles of surrogates in mesh adaptive direct search
    *Journal of Global Optimization*, 70(3):645–675, 2018.

  .. [AuLedSa2021] C.Audet, S.Le Digabel and R.Saltet.
    Quantifying uncertainty with ensembles of surrogates for blackbox optimization.
    Rapport technique G-2021-37, Les cahiers du GERAD, 2021.
    http://www.optimization-online.org/DB_HTML/2021/07/8489.html