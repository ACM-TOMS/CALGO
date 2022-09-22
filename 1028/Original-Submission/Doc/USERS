----------------------------------
VTMOP Organization and Usage Guide
----------------------------------

VTMOP offers two interfaces for solving blackbox MOPs. They are
 * the return-to-caller interface and
 * the driver subroutine VTMOP_SOLVE.

In many situations, the objective function F can easily be wrapped
in a Fortran subroutine, and the driver subroutine interface is preferred
because of its ease of use. However, the return-to-caller interface offers
a flexible alternative for situations where the computing environment
makes wrapping F in a Fortran subroutine difficult or inefficient.
Both interfaces are implemented in ISO Fortran 2008 and support some forms
of OpenMP parallelism.

For information on how to implement your objective function as a Fortran
function, including cases where F is an ISO C/C++ function and when
F is a command line executable, see OBJ_FUNC_README. 

Further usage information is provided in the comments at the top of
each subroutine, in their respective source files. A list of source files
and their contents is provided in the Organization section of this
document.

The VTMOP_SOLVE Driver
----------------------

The driver subroutine VTMOP_SOLVE implements the same process as in
the Return-to-Caller interface (below) but without reading and writing
the state and using a Fortran implementation of F.
For more information on how VTMOP performs these steps, see the TOMS Algorithm
paper.

The search phase is performed either via DIRECT or via Latin hypercube
design, with the choice being specified as an optional input.
For a search via DIRECT, the search budget inputs are used to specify
the iteration limit for VTDIRECT95. For a search via Latin hypercube
design, the search budgets are used to specify the number of points in
the design. In addition, VTMOP_SOLVE offers a checkpointing system and
a parallel option.

Remark U.1
 * When the return-to-caller interface is used, the user is expected to
   handle all function evaluation data. Indeed, this is the motivation for
   using the return-to-caller interface. Therefore, the user is also
   responsible for saving function evaluation data for checkpointing.
   The checkpointing system in VTMOP separately saves iteration
   data (such as the sequence of local trust regions (LTRs) and parameter
   settings) and function evaluation data so that the iteration data can
   still be tracked when the return-to-caller interface is used, although
   the function evaluation data is tracked only when using VTMOP_SOLVE.
   As a side effect, even when using VTMOP_SOLVE, the checkpoint files do
   not track which design point evaluations have been requested but not
   fulfilled. So, if an instance of VTMOP_SOLVE terminates midway through
   a batch of function evaluations, reloading the last checkpoint may
   result in a new sequence of evaluations.

The Return-to-Caller Interface
------------------------------

In many real-world blackbox optimization problems, special
purpose computing environments and libraries for coordinating the use
of parallel resources require the evaluation of F to be decoupled
from the optimization algorithm. In these situations, it is convenient
to offer a return-to-caller interface, where F is evaluated only
from outside of any Fortran subroutine's call stack.
Here it is necessary to preserve the entire state of VTMOP's databases
and internal variables in between calls to F. The state is recorded
in a Fortran derived data type VTMOP_TYPE, created by the subroutine
VTMOP_INIT, and stored between invocations of components of VTMOP.
The cycle is as follows:
 1. k -> 0;
 2. CALL VTMOP_INIT (this creates a VTMOP_TYPE structure holding the problem
                     parameters and VTMOP's state variables);
 3. CALL VTMOP_LTR (this computes the lower and upper bounds for the next LTR);
 4. evaluate F at user-selected x-values, which are space-filling in
    the LTR returned by VTMOP_LTR;
 5. CALL VTMOP_OPT (this computes the kth batch of candidate designs);
 6. evaluate the kth batch of candidate designs from VTMOP_OPT with F;
 7. k -> k+1;
 8. check the termination conditions (max budget or iteration limit),
    and either terminate by calling VTMOP_FINALIZE or return to Step 3.

Parallel Interface
------------------

Two opportunities for parallelism are offered in VTMOP.
 * First, the iteration tasks in the worker subroutines VTMOP_LTR and
   VTMOP_OPT can be parallelized. This includes computing the Delaunay
   graph in parallel using the parallel DELAUNAYSPARSEP driver in
   delsparse.f90, fitting the surrogates in parallel, and solving the
   all surrogate problems in parallel. This iteration level parallelism
   is provided via OpenMP 4.5 and is available through either the
   return-to-caller interface or the driver subroutine.
 * The second source of parallelism is concurrently evaluating F at points
   requested by VTMOP. In the context of blackbox optimization, this is
   typically the more impactful source of parallelism. Since the
   return-to-caller interface does not evaluate F, this form of parallelism
   is available only through the driver subroutine (which offers a choice
   between iteration task parallelism, concurrent function evaluations,
   or both). Recall that VTMOP requires function evaluations only while
   searching LTRs and when evaluating candidate designs. In both cases,
   concurrent evaluations are achieved by spawning a new OpenMP task for
   each evaluation of F. Note that because the actual function evaluations
   are handled by the user, the user is also responsible for capping the
   number of actual function evaluations in order to prevent the
   oversubscription of resources.

Remark U.2
 * Within a single evaluation of F, there may be ample opportunities for
   parallelism, which is problem dependent and left to the user.
   In particular, in many real-world applications, each evaluation of F
   could depend on output from a multinode, distributed simulation.
   VTMOP does not distribute calls to F, but it does provide the
   ability for concurrent evaluations of F via OpenMP task-based parallelism.
   This places the burden of distributing calls to F on the user, but it
   allows for greater flexibility.

Remark U.3
 * When evaluating candidate designs and when using the Latin hypercube
   design during the search phase, all candidates/design points can be
   evaluated in parallel if enough computing resources are available.
   For a search via DIRECT, VTMOP uses a slightly modified implementation
   of VTdirect, called bVTdirect, which uses OpenMP tasks to perform
   concurrent evaluations of box centers. During a search phase, VTMOP makes
   several asynchronous calls to bVTdirect, resulting in nested parallelism.

Remark U.4
 * The following considerations apply only when using the DIRECT based search
   in the driver VTMOP_SOLVE: As described above, running asynchronous
   instances of DIRECT will result in requests for redundant design point
   evaluations, which VTMOP handles by using OpenMP task dependency.
   This means that when using bVTdirect during the search phase, many
   OpenMP tasks will be waiting on the completion of other tasks.
   As a result, the total number of OpenMP threads should be greater than
   the desired number of concurrent function evaluations so that tasks that
   are waiting will not block other designs from being evaluated.
   A general recommendation is to set the number of OpenMP level one threads
   to the desired number of concurrent function evaluations.
   (For example, if each evaluation of F requires a single node, then the
   number of level one threads should be the number of available nodes).
   This number will also be used for determining the number of concurrent
   evaluations when evaluating the candidate designs and for determining the
   number of iteration tasks that can be performed in parallel.
   If the number of total resources is greater than the number of objectives,
   then the number of level-two threads should typically be two or three.
   If the number of objectives is greater than the amount of resources,
   then a larger number of level two threads may be needed.
   In OpenMP 4.5, these values can be controlled using the OMP_NUM_THREADS
   environment variable. Use the bash command:
   
   export OMP_NUM_THREAD=[number of lvl 1 threads],[number of lvl 2 threads]

   In order to allow for nested parallelism, which is needed so that
   multiple instances of bVTdirect can run at once, the environment variable
   OMP_NESTED should be set to true. Use the bash command:
   
   export OMP_NESTED=TRUE

Remark U.5
 * Each serial run of VTMOP_SOLVE is deterministic when using
   the search via DIRECT, LSHEP surrogate, and generalized pattern search
   (GPS) optimization options. For grid aligned data, however, LSHEP can
   be sensitive to the ordering of function evaluations in VTMOP's
   internal databases. Since DIRECT samples on an implicit mesh, using
   concurrent function evaluations in a run of VTMOP_SOLVE with the
   above settings could produce nondeterministic (but still reasonable)
   results.

Organization
------------

The package VTMOP is organized as follows.
 * The file vtmop.f90 contains the subroutines VTMOP_SOLVE and 
   DELAUNAYGRAPH and the modules VTMOP_MOD and VTMOP_LIB.
   VTMOP_MOD contains all of the worker subroutines described in the TOMS
   Algorithm paper and several others that are used by VTMOP_SOLVE.
   VTMOP_LIB provides interfaces for VTMOP_SOLVE, DELAUNAYGRAPH, and
   those subroutines from VTMOP_MOD that are externally useful,
   including VTMOP_INIT, VTMOP_LTR, VTMOP_OPT, and VTMOP_FINALIZE.
 * The file vtmop_func.f90 provides implementations of several analytic
   multiobjective test functions from the literature, which match the
   interface to the user written subroutine OBJ_FUNC, defining F.
 * The file sample.f90 provides sample code for minimizing one of the
   multiobjective test problems in vtmop_func.f90 using VTMOP_SOLVE
   and verifying the output.
 * The file sVTdirect.f90 contains a lightweight version of the serial
   driver from VTDIRECT95, with features that are irrelevant to VTMOP removed.
 * The file bVTdirect.f90 contains the modification of sVTdirect.f90 described
   above.
 * The file shared_modules.f90 contains the shared modules required
   by {s|b}VTdirect.f90. Notably, the module REAL_PRECISION provides
   access to the real kind R8, which is used for approximately 64-bit
   precision on all known machines.
 * The files linear_shepard.f90, qnstop.f90, and delsparse.f90 respectively
   provide subsets of the packages SHEPPACK, QNSTOP, and DELAUNAYSPARSE
   that are relevant to VTMOP.
 * Additionally, the files blas.f, lapack.f, and slatec.f respectively
   provide subsets of BLAS, LAPACK, and SLATEC that are used by VTMOP
   or its dependencies. If any of these libraries are already installed,
   link against those copies instead.

The dependency chart for these files is further detailed in depend_graph.txt.

Remark U.6
 * The worker subroutine VTMOP_INIT accepts numerous inputs
   (many optional) when initializing the data type VTMOP_TYPE.
   Most of these inputs are further discussed in the TOMS Algorithm paper,
   and they are summarized here. The required inputs are D and P,
   which specify the number of design variables and objectives, respectively;
   and LB and UB, which specify the lower and upper bound constraints,
   respectively.

   The optional inputs are:
    * DES_TOL and OBJ_TOL, which specify the design space tolerance and
      objective space tolerance, respectively (by default, the square root
      of the working precision EPS, defined below);
    * TRUST_RADF, DECAY, and MIN_RADF, which specify the initial trust
      region radius (as a fraction of UB-LB), the decay rate, and the
      minimum trust region radius (as a fraction of UB-LB), respectively
      (by default, 20% of each dimension of the feasible design space,
      0.5, and 10% of the initial trust region radius fraction); 
    * EPS, which specifies the working precision of the machine
      (by default, the square root of the unit round-off);
    * EPSW, which specifies the fudge factor for zero weights
      (by default, the fourth root of the unit round-off);
    * LOCAL_OPT and LOPT_BUDGET, which specify the procedure for performing
      local optimization and its iteration budget, respectively
      (by default, GPS with a budget of 2,500 iterations);
    * FIT_SURROGATES and EVAL_SURROGATES, which specify the procedures
      for fitting and evaluating the surrogates, respectively
      (by default, wrappers for the subroutines LSHEP and LSHEP_VAL
      from SHEPPACK);
    * OBJ_BOUNDS, which specifies lower and upper bounds on the interesting
      objective ranges (by default, there are no bounds);
    * PFLAG, which is a Boolean value specifying whether to perform parallel
      iteration tasks (by default, false); and
    * ICHKPT, which specifies the checkpointing mode (i.e., no checkpointing,
      use checkpointing, or recover from checkpoint, with no checkpointing
      by default). When checkpointing is used, the return-to-caller interface
      stores its iteration data in vtmop.chkpt, and the user is responsible
      for saving function evaluation data.

   Every VTMOP subroutine also returns an integer valued error flag IERR,
   which contains zero after a successful operation or a three digit error
   code if a failure occurs. The specific error codes are contained in
   the comments above each subroutine definition.

Remark U.7
 * The driver VTMOP_SOLVE accepts all of the inputs described in Remark
   U.6 except PFLAG, which is replaced by an integer PMODE, specifying
   whether to use iteration parallelism, parallelism between function
   evaluations, or both (by default, neither). Also, when checkpointing
   is used, VTMOP_SOLVE stores function evaluation data in a second
   checkpointing file, vtmop.dat.
   In addition to the inputs in Remark U.6, VTMOP_SOLVE requires a
   subroutine OBJ_FUNC that evaluates the blackbox multiobjective
   function F; and VTMOP_SOLVE returns EFFICIENT_X and
   PARETO_F, containing the efficient and nondominated sets, respectively.
   Additional optional arguments for VTMOP_SOLVE are as follows.
    * ADAPTIVE_SEARCH is a Boolean variable that specifies whether the
      search via DIRECT should be used when true, as opposed to the
      Latin hypercube search when false (by default, true).
    * BB_BUDGET is the limit on the number of blackbox function evaluations
      (i.e., evaluations of F, by default, 1,000).
    * MAXITERS is the iteration limit for VTMOP_SOLVE (by default, unlimited).
    * INITIAL_SBUDGET and SEARCH_BUDGET are the budget per search
      when k=0 and k>0, respectively -- when ADAPTIVE_SEARCH is true,
      the search budgets specify the number of iterations for DIRECT
      (by default, 10 and 5, respectively); otherwise, they specify the
      number of points in each Latin hypercube design (by default, 16*D^2
      and 8*D, respectively).
    * On input, DES_PTS and OBJ_PTS can be used to specify the
      initial databases of precomputed function values; and on output,
      they will be reallocated and used to return the final database
      of all function evaluations.

Remark U.8
 * VTMOP allows for certain optional inputs to be changed after recovering
   from a checkpoint. Therefore, the following optional inputs are not
   saved when checkpointing is active; they must be resupplied by the user
   if nondefault values are required:
    * LOPT_BUDGET,
    * LOCAL_OPT,
    * FIT_SURROGATES,
    * EVAL_SURROGATES,
    * PMODE/PFLAG,
    * ADAPTIVE_SEARCH,
    * BB_BUDGET,
    * MAXITERS,
    * INITIALS_BUDGET,
    * SEARCH_BUDGET,
    * DES_PTS, and
    * OBJ_PTS.

Remark U.9
 * The cost of each search phase can get out of hand if one is not
   mindful of the interplay between INITIAL_SBUDGET and SEARCH_BUDGET and
   D when using the setting ADAPTIVE_SEARCH=.TRUE.
   For large values of D, setting ADAPTIVE_SEARCH=.FALSE. is the better choice.
