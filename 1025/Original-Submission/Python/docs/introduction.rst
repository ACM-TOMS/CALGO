Introduction
============

We consider a general minimization problem:

.. math::
   \min_\mathbf{x} \, y(\mathbf{x})



Bayesian optimization proceeds through construction of a surrogate cost
function :math:`\tilde{y}(\mathbf{x})`. This surrogate is represented as a basis
function expansion, around each evaluated point(:math:`\mathbf{x}_i, i \in [1,N]`).
This ensures that the surrogate passes through (interpolates) the evaluated
points. In the case of evaluations with noisy data, the surrogate shall pass
within 1 standard deviation from the mean at the evaluated points.
Analytically, the surrogate :math:`\tilde{y}(\mathbf{x})` after N function
evaluations is represented as

.. math::
    \tilde{y}(\mathbf{x}) = \sum_{i=1}^N w_i k(\mathbf{x}, \mathbf{x}_i)

where :math:`k(\mathbf{x}, \mathbf{x}_i)` is a kernel function, i.e., it takes in two arguments,
:math:`\mathbf{x}, \, \mathbf{x}_i` and returns a scalar value. This scalar is representative of how correlated is
the function :math:`y(\mathbf{x})` at :math:`\mathbf{x}` and :math:`\mathbf{x}_i`. The weights :math:`w_i` are
calculated by solving the system of :math:`N` linear equations in :math:`w_i`. In matrix notation, this is represented
using a covariance matrix(:math:`\mathbf{K}`):

.. math::
    \mathbf{K}\, \bar{w} & = y \\
    \mathbf{K}_{i,j} & = k(\mathbf{x}_i, \mathbf{x}_j) , \, \, i,j\in[1,N] \\
    y_i & = y(\mathbf{x}_i) , \, \,i\in[1,N] \\
    \bar{w} & = \{w_i\}, \, \,i\in[1,N]


Hence the weights are calculated through the inversion :math:`\bar{w} = \mathbf{K}^{-1}\,y`. Note that the covariance
matrix :math:`\mathbf{K}` is a Gram matrix of a positive definite kernel function, making it symmetric and positive
semi-definite. Furthermore, since with every iteration only a finite number of rows are added to the covariance matrix,
efficient inversion is possible through incremental Cholesky decomposition. The mean and variance of the surrogate
are then calculated as:

.. math::
    \mu(\mathbf{x}_{N+1}) & = \mathbf{k}^T \mathbf{K}^{-1} y_{1:N} \\
    \sigma^2(\mathbf{x}_{N+1}) & = k(\mathbf{x}_{N+1}, \mathbf{x}_{N+1}) - \mathbf{k}^T\,\mathbf{K}^{-1}\,\mathbf{k}

where

.. math::
    \mathbf{k} = k(\mathbf{x}_{1:N}, \mathbf{x}_{N+1}) = [k(\mathbf{x}_1,\mathbf{x}_{N+1})\, k(\mathbf{x}_2,\mathbf{x}_{N+1})\, . . .            k(\mathbf{x}_N,\mathbf{x}_{N+1})]

At each iteration, the surrogate is updated with new data from the cost function. The locations where the next
evaluation is done is determined through optimization of an *acquisition function*. An acquisition function is a means
to estimate the new information content at a location. It uses the mean and variance calculated in the above steps.

Some sample radial kernel functions include:

* squared exponential kernel function : Infinitely differentiable

.. math::
   k(r) = \theta_0 exp\Bigg(- \frac{r^2}{\theta^2}\Bigg)

* Matern class of kernel function :

.. math::
   k_{Matern}(r) = \frac{2^{1-\nu}}{\Gamma(\nu)}\Bigg(\frac{\sqrt{2\nu}r}{l}\Bigg)^\nu K_\nu\Bigg(\frac{\sqrt{2\nu}r}{l}\Bigg)

where :math:`K_\nu` is the modified Bessel function, :math:`\nu,l` are positive constants

* Rational quadratic kernel function:

.. math::
   k_{RQ}(r) = \Bigg(1 + \frac{r^2}{2\alpha l^2}\Bigg)^{-\alpha}


where :math:`r = ||\mathbf{x}_1 - \mathbf{x}_2||`


Some example acquisition functions are:

* Confidence bounds

.. math::
   LCB = \mu - \kappa \sigma

* Probability of improvement

.. math::
   PI = \mathbf{cdf}(\gamma)


* Expectation of improvement

.. math::
   EI = sqrt(variance) * (\gamma * \mathbf{cdf}(\gamma) + \mathbf{pdf}(\gamma))

where :math:`\gamma = \frac{\mu}{\sigma}`, :math:`\mathbf{cdf}` is cumulative normal distribution function and
:math:`\mathbf{pdf}` is normal probability distribution function
